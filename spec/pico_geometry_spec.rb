# frozen_string_literal: true

require "spec_helper"

# Tier 1: pure 2D geometry tests for the Pico inferior-2/3 method.
#
# These exercise GlenoidMorphology::Measure helpers directly on synthetic 2D
# rim point clouds. No 3D mask, no NIfTI affine -- this isolates the Pico
# fit and the inferior-2/3 filter from the rest of the pipeline.
RSpec.describe "GlenoidMorphology::Measure (Pico 2D geometry)" do
  M = GlenoidMorphology::Measure

  # Generate rim points on a circle of radius `r` over an angular range,
  # measured CCW from the +x axis.
  def arc_points(r:, center: [0.0, 0.0], from_deg:, to_deg:, n:)
    pts = Numo::DFloat.zeros(n, 2)
    n.times do |i|
      t = from_deg + (to_deg - from_deg) * i.to_f / [n - 1, 1].max
      rad = t * Math::PI / 180.0
      pts[i, 0] = center[0] + r * Math.cos(rad)
      pts[i, 1] = center[1] + r * Math.sin(rad)
    end
    pts
  end

  def concat(*arrs)
    total = arrs.sum { |a| a.shape[0] }
    out = Numo::DFloat.zeros(total, 2)
    off = 0
    arrs.each do |a|
      n = a.shape[0]
      out[off...(off + n), true] = a
      off += n
    end
    out
  end

  describe ".inferior_two_thirds_mask" do
    it "keeps points within +/- 120° of the inferior pole" do
      # Inferior pole is -y in this synthetic frame.
      si_dir = Numo::DFloat[0.0, -1.0]
      centroid = Numo::DFloat[0.0, 0.0]
      # 36 evenly-spaced rim points, every 10°.
      pts = arc_points(r: 10.0, from_deg: 0.0, to_deg: 350.0, n: 36)
      mask = M.inferior_two_thirds_mask(pts, centroid: centroid, si_dir2d: si_dir)
      # Expect ~24 of 36 points (240° / 360°).
      expect(mask.count_true).to be_within(2).of(24)
    end

    it "puts the inferior-pole point inside the mask" do
      si_dir = Numo::DFloat[0.0, -1.0]
      centroid = Numo::DFloat[0.0, 0.0]
      pts = Numo::DFloat[[0.0, -10.0], [0.0, 10.0]]
      mask = M.inferior_two_thirds_mask(pts, centroid: centroid, si_dir2d: si_dir)
      expect(mask[0]).to eq(1) # inferior pole IN
      expect(mask[1]).to eq(0) # superior pole OUT
    end
  end

  describe ".longest_pca_axis_2d" do
    it "recovers an axis-aligned long axis" do
      # A vertical bar of points: longest axis should be (0, ±1).
      n = 50
      pts = Numo::DFloat.zeros(n, 2)
      n.times { |i| pts[i, 1] = (i - n / 2) * 1.0; pts[i, 0] = (i.odd? ? 0.1 : -0.1) }
      axis = M.longest_pca_axis_2d(pts)
      expect(axis[0].abs).to be < 0.1
      expect(axis[1].abs).to be_within(0.05).of(1.0)
    end
  end

  describe ".fit_rim_circle (on the Pico subset)" do
    it "recovers ~12.5mm radius when the lower arc is clean and the upper arc is dented inward by 4mm" do
      # Inferior 2/3 = clean circular arc, r = 12.5, theta in [-120°, 120°] of inferior pole.
      # We'll put inferior pole at -y, so the clean arc is the lower portion.
      # In CCW angles (from +x axis): inferior pole = -90°.
      # |theta_from_inferior| <= 120° <=> angle in CCW from +x in [-210°, 30°],
      # i.e. wrap to [150°..360°] U [0°..30°].
      # Easier: build the inferior arc as theta in [-90 - 120, -90 + 120] = [-210, 30]°.
      lower = arc_points(r: 12.5, from_deg: -210.0, to_deg: 30.0, n: 240)
      # Superior 1/3 = inward-dented arc: theta in [30°, 150°], but at r=8.5 (4mm inward).
      upper = arc_points(r: 8.5, from_deg: 30.0, to_deg: 150.0, n: 120)
      all_pts = concat(lower, upper)

      opts = GlenoidMorphology::Measure::DEFAULT_OPTS.merge(seed: 42)
      si_dir = Numo::DFloat[0.0, -1.0]
      centroid = all_pts.mean(axis: 0)
      mask = M.inferior_two_thirds_mask(all_pts, centroid: centroid, si_dir2d: si_dir)
      pico_pts = all_pts[mask.where, true]

      circle = M.fit_rim_circle(pico_pts, opts)
      expect(circle.radius).to be_within(0.5).of(12.5)
    end

    it "recovers ~12.5mm radius when the inferior arc has a 30° chord-gap" do
      # Inferior 2/3 but with a 30° gap centred at the inferior pole.
      # Inferior pole CCW = -90°. So we have arcs [-210, -105] and [-75, 30].
      seg_a = arc_points(r: 12.5, from_deg: -210.0, to_deg: -105.0, n: 110)
      seg_b = arc_points(r: 12.5, from_deg: -75.0,  to_deg: 30.0,   n: 110)
      upper = arc_points(r: 8.5,  from_deg: 30.0,   to_deg: 150.0,  n: 120)
      all_pts = concat(seg_a, seg_b, upper)

      opts = GlenoidMorphology::Measure::DEFAULT_OPTS.merge(seed: 42)
      si_dir = Numo::DFloat[0.0, -1.0]
      centroid = all_pts.mean(axis: 0)
      mask = M.inferior_two_thirds_mask(all_pts, centroid: centroid, si_dir2d: si_dir)
      pico_pts = all_pts[mask.where, true]
      circle = M.fit_rim_circle(pico_pts, opts)
      expect(circle.radius).to be_within(0.6).of(12.5)
      # Centre should be near origin (the rim was drawn around origin).
      expect(circle.cx.abs).to be < 1.0
      expect(circle.cy.abs).to be < 1.0
    end
  end
end

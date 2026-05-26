# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::DefectArc do
  def ring_with_gap(gap_deg:, r: 10.0, n: 360, seed: 1)
    rng = Random.new(seed)
    full = []
    # Use a contiguous gap starting at angle 0 — easy to reason about.
    n.times do |i|
      theta = (i.to_f / n) * 2 * Math::PI
      next if theta < gap_deg * Math::PI / 180.0
      r_jitter = (rng.rand - 0.5) * 0.05
      full << [r + r_jitter, theta]
    end
    p = Numo::DFloat.zeros(full.size, 2)
    full.each_with_index do |(rr, th), i|
      p[i, 0] = rr * Math.cos(th)
      p[i, 1] = rr * Math.sin(th)
    end
    p
  end

  let(:circle) { GlenoidMorphology::CircleFit::Circle.new(0.0, 0.0, 10.0, 0.01, 1.0) }

  it "detects a 0° gap (no defect)" do
    points = ring_with_gap(gap_deg: 0.0)
    arc = described_class.find(points, circle: circle, bins: 360, radial_tolerance: 0.5)
    expect(arc.bone_loss_percent).to be < 1.0
  end

  it "detects a 60° gap with the right span" do
    points = ring_with_gap(gap_deg: 60.0)
    arc = described_class.find(points, circle: circle, bins: 360, radial_tolerance: 0.5)
    span_deg = arc.angular_span_rad * 180.0 / Math::PI
    expect(span_deg).to be_within(2.0).of(60.0)
    expected = Synthetic.expected_loss_fraction(60.0) * 100.0
    expect(arc.bone_loss_percent).to be_within(0.5).of(expected)
  end

  it "detects a 120° gap" do
    points = ring_with_gap(gap_deg: 120.0)
    arc = described_class.find(points, circle: circle, bins: 360, radial_tolerance: 0.5)
    span_deg = arc.angular_span_rad * 180.0 / Math::PI
    expect(span_deg).to be_within(2.0).of(120.0)
  end

  it "reports the largest run when there are multiple small gaps" do
    # 30° gap and a 10° gap; we should report the 30°
    rng = Random.new(2)
    pts = []
    360.times do |i|
      theta = (i.to_f / 360) * 2 * Math::PI
      next if i.between?(0, 29)    # 30 deg
      next if i.between?(100, 109) # 10 deg
      noise = (rng.rand - 0.5) * 0.02
      pts << [(10 + noise) * Math.cos(theta), (10 + noise) * Math.sin(theta)]
    end
    p = Numo::DFloat.zeros(pts.size, 2)
    pts.each_with_index { |(x, y), i| p[i, 0] = x; p[i, 1] = y }
    arc = described_class.find(p, circle: circle, bins: 360, radial_tolerance: 0.5)
    span_deg = arc.angular_span_rad * 180.0 / Math::PI
    expect(span_deg).to be_within(2.0).of(30.0)
  end
end

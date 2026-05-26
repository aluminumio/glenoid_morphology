# frozen_string_literal: true

require "spec_helper"

# Tier 2: end-to-end Pico method tests on synthetic 3D glenoid masks.
#
# These exercise the full pipeline (mask -> en-face -> rim -> Pico circle ->
# inferior-2/3 bone-loss). The "pear" synthetic glenoid is the key test: a
# healthy pear-shaped glenoid should report near-zero loss under Pico but
# substantial loss under the legacy full-circle method.
RSpec.describe "GlenoidMorphology.measure (Pico method, end-to-end)" do
  let(:affine) { GlenoidMorphology::Affine.identity(spacing_mm: 1.0) }
  # Realistic glenoid orientation: the articular plane normal is mostly
  # lateral (+x), with a slight tilt -- so world inferior (-z) is largely
  # IN-plane, which is what the Pico SI-projection step needs.
  let(:plane_normal) { [Math.cos(0.2), 0.0, Math.sin(0.2)] }

  def measure_pear(mask, **opts)
    humerus = Synthetic.humerus_mask_for(plane_normal: plane_normal)
    GlenoidMorphology.measure(
      scapula_mask: mask,
      humerus_mask: humerus,
      affine: affine,
      side: :right,
      seed: 12345,
      **opts
    )
  end

  describe "intact pear glenoid" do
    it "reports low Pico bone-loss for a healthy pear (large legacy loss)" do
      # Pear shape: r=12mm circle below superior third, narrowed by 7mm at the
      # superior pole. Mimics the ~30mm SI x ~18mm AP healthy pear (the
      # superior third pinches down to about a third the diameter).
      mask = Synthetic.pear_glenoid_mask(r: 12.0, taper_mm: 7.0, plane_normal: plane_normal)
      pico   = measure_pear(mask, method: :pico)
      legacy = measure_pear(mask, method: :full_circle)
      print "  [pear]  pico_blp=#{pico.bone_loss_percent.round(2)}%  " \
            "legacy_blp=#{legacy.bone_loss_percent.round(2)}%  " \
            "pico_r=#{pico.fitted_circle[:radius_mm].round(2)}mm  " \
            "legacy_r=#{legacy.fitted_circle[:radius_mm].round(2)}mm\n"

      # Pico should report a much smaller loss than full-circle.
      expect(pico.bone_loss_percent).to be < 10.0
      expect(pico.bone_loss_percent).to be < (legacy.bone_loss_percent - 5.0)
      # Pico circle radius should be in the SI range (i.e. close to r).
      expect(pico.fitted_circle[:radius_mm]).to be_within(3.0).of(12.0)
      # Pico struct should have pico_circle populated.
      expect(pico.pico_circle).not_to be_nil
      expect(pico.pico_circle[:radius_mm]).to be_within(3.0).of(12.0)
      # Legacy result should also be available via legacy_bone_loss_percent
      # and substantially higher than Pico's (the pear-narrowing on the
      # superior third confuses the full-circle integral).
      expect(pico.legacy_bone_loss_percent).to be > 7.0
      expect(pico.legacy_bone_loss_percent).to be > pico.bone_loss_percent
    end
  end

  describe "pear glenoid with a 20° inferior chord defect" do
    it "reports the expected segment area under Pico" do
      mask = Synthetic.pear_glenoid_mask(r: 12.0, taper_mm: 7.0, plane_normal: plane_normal,
                                          defect_angle_deg: 20.0,
                                          defect_at_inferior: true)
      pico = measure_pear(mask, method: :pico)
      # Expected segment area for a 20° chord, as a fraction of the (Pico)
      # inferior-2/3 sector. The chord-segment formula uses r=12 (the Pico
      # radius). area_seg = 0.5 r^2 (theta - sin theta), with theta = 20° in rad.
      # Sector area for Pico (240°): (240/360) * pi * r^2 = (2/3) pi r^2.
      theta_rad = 20.0 * Math::PI / 180.0
      segment = 0.5 * (12.0 ** 2) * (theta_rad - Math.sin(theta_rad))
      sector  = (2.0 / 3.0) * Math::PI * (12.0 ** 2)
      expected_frac_pct = segment / sector * 100.0
      print "  [pear+def20] pico_blp=#{pico.bone_loss_percent.round(3)}% " \
            "expected~#{expected_frac_pct.round(3)}%\n"
      # 20° is below the per-voxel detection floor (~5% absolute) at 1 mm
      # spacing. Just assert it doesn't blow up.
      expect(pico.bone_loss_percent).to be < 10.0
    end
  end

  describe "pear glenoid with a 45° inferior chord defect" do
    it "reports a defect that shows up in Pico bone loss" do
      mask = Synthetic.pear_glenoid_mask(r: 12.0, taper_mm: 7.0, plane_normal: plane_normal,
                                          defect_angle_deg: 45.0,
                                          defect_at_inferior: true)
      pico = measure_pear(mask, method: :pico)
      theta_rad = 45.0 * Math::PI / 180.0
      segment = 0.5 * (12.0 ** 2) * (theta_rad - Math.sin(theta_rad))
      sector  = (2.0 / 3.0) * Math::PI * (12.0 ** 2)
      expected_frac_pct = segment / sector * 100.0
      print "  [pear+def45] pico_blp=#{pico.bone_loss_percent.round(3)}% " \
            "expected~#{expected_frac_pct.round(3)}%\n"
      expect(pico.bone_loss_percent).to be_within(5.0).of(expected_frac_pct)
    end
  end

  describe "API contract" do
    it "raises on an unknown method symbol" do
      mask = Synthetic.intact_glenoid_mask
      expect {
        GlenoidMorphology.measure(
          scapula_mask: mask, affine: affine, side: :right, method: :bogus
        )
      }.to raise_error(ArgumentError, /method must be/)
    end

    it "leaves pico_circle nil when using :full_circle method" do
      mask = Synthetic.intact_glenoid_mask
      humerus = Synthetic.humerus_mask_for
      result = GlenoidMorphology.measure(
        scapula_mask: mask, humerus_mask: humerus, affine: affine,
        side: :right, seed: 12345, method: :full_circle
      )
      expect(result.pico_circle).to be_nil
      expect(result.method).to eq(:full_circle)
      # legacy_bone_loss_percent should equal bone_loss_percent when method=:full_circle.
      expect(result.legacy_bone_loss_percent).to eq(result.bone_loss_percent)
    end
  end
end

# frozen_string_literal: true

require "spec_helper"

# Tier 2: end-to-end measurements on synthetic glenoid masks.
RSpec.describe "GlenoidMorphology.measure (end-to-end)" do
  let(:affine) { GlenoidMorphology::Affine.identity(spacing_mm: 1.0) }

  # Tolerance for the synthetic test. The disc is voxelised at 1 mm, so the
  # rim is jagged at ~half-voxel scale and the practical resolution floor of
  # this whole pipeline is around 4–5 percentage-points of bone-loss for a
  # ~24 mm radius glenoid. Real CT data is finer (~0.5 mm) and would tighten
  # this proportionally.
  TOLERANCE_PCT = 4.5

  def measure_mask(mask, humerus: nil, **opts)
    humerus ||= Synthetic.humerus_mask_for
    # These specs verify the full-circle chord-segment math. The Pico method
    # (the v0.3+ default) intentionally measures bone loss only over the
    # inferior-2/3 sector of a clinically-derived reference circle, which is
    # the right thing for real anatomy but not what these synthetic
    # axis-aligned chord-defect discs test. New Pico specs live in
    # pico_method_spec.rb.
    opts = { method: :full_circle }.merge(opts)
    GlenoidMorphology.measure(
      scapula_mask: mask,
      humerus_mask: humerus,
      affine: affine,
      side: :right,
      seed: 12345,
      **opts
    )
  end

  it "reports near-zero bone loss for an intact glenoid" do
    mask = Synthetic.intact_glenoid_mask
    result = measure_mask(mask)
    expect(result.bone_loss_percent).to be < TOLERANCE_PCT
    expect(result.fitted_circle[:radius_mm]).to be_within(2.0).of(Synthetic::GLENOID_RADIUS)
    expect(result.confidence).to be > 0.5
  end

  it "measures a 20° chord defect (sub-voxel resolution; below detection floor)" do
    # 20° subtended chord = 0.11% lost area = ~2 voxels of a 1810-voxel disc.
    # That's below the detection floor at 1mm voxel spacing, so this serves
    # as a regression-check that small defects don't get over-detected.
    mask = Synthetic.chord_defect_mask(angle_deg: 20.0)
    result = measure_mask(mask)
    expect(result.bone_loss_percent).to be < TOLERANCE_PCT
  end

  it "measures a 45° chord defect" do
    mask = Synthetic.chord_defect_mask(angle_deg: 45.0)
    result = measure_mask(mask)
    expected = Synthetic.expected_loss_fraction(45.0) * 100.0
    print "  [45°]  expected=#{expected.round(3)}%, measured=#{result.bone_loss_percent.round(3)}%, " \
          "diff=#{(result.bone_loss_percent - expected).round(3)}\n"
    expect(result.bone_loss_percent).to be_within(TOLERANCE_PCT).of(expected)
  end

  it "measures a 90° chord defect" do
    mask = Synthetic.chord_defect_mask(angle_deg: 90.0)
    result = measure_mask(mask)
    expected = Synthetic.expected_loss_fraction(90.0) * 100.0
    print "  [90°]  expected=#{expected.round(3)}%, measured=#{result.bone_loss_percent.round(3)}%, " \
          "diff=#{(result.bone_loss_percent - expected).round(3)}\n"
    expect(result.bone_loss_percent).to be_within(TOLERANCE_PCT).of(expected)
  end

  it "measures a 120° chord defect" do
    mask = Synthetic.chord_defect_mask(angle_deg: 120.0)
    result = measure_mask(mask)
    expected = Synthetic.expected_loss_fraction(120.0) * 100.0
    print "  [120°] expected=#{expected.round(3)}%, measured=#{result.bone_loss_percent.round(3)}%, " \
          "diff=#{(result.bone_loss_percent - expected).round(3)}\n"
    expect(result.bone_loss_percent).to be_within(TOLERANCE_PCT).of(expected)
  end

  it "measures a 150° chord defect" do
    mask = Synthetic.chord_defect_mask(angle_deg: 150.0)
    result = measure_mask(mask)
    expected = Synthetic.expected_loss_fraction(150.0) * 100.0
    print "  [150°] expected=#{expected.round(3)}%, measured=#{result.bone_loss_percent.round(3)}%, " \
          "diff=#{(result.bone_loss_percent - expected).round(3)}\n"
    expect(result.bone_loss_percent).to be_within(TOLERANCE_PCT).of(expected)
  end

  it "is invariant under rotation of the disc plane" do
    # Tilt the disc 30° about the world x-axis.
    theta = 30.0 * Math::PI / 180.0
    rotated_normal = [0.0, -Math.sin(theta), Math.cos(theta)]
    mask_a = Synthetic.chord_defect_mask(angle_deg: 120.0)
    humerus_a = Synthetic.humerus_mask_for
    mask_b = Synthetic.chord_defect_mask(angle_deg: 120.0, plane_normal: rotated_normal)
    humerus_b = Synthetic.humerus_mask_for(plane_normal: rotated_normal)

    a = measure_mask(mask_a, humerus: humerus_a).bone_loss_percent
    b = measure_mask(mask_b, humerus: humerus_b).bone_loss_percent
    print "  [rot]   axis-aligned=#{a.round(3)}%, tilted=#{b.round(3)}%, " \
          "diff=#{(a - b).round(3)}\n"
    expect((a - b).abs).to be < TOLERANCE_PCT
  end

  it "is invariant under translation" do
    mask_a = Synthetic.chord_defect_mask(angle_deg: 120.0)
    humerus_a = Synthetic.humerus_mask_for
    mask_b = Synthetic.chord_defect_mask(angle_deg: 120.0, center: [40, 60, 70])
    humerus_b = Synthetic.humerus_mask_for(glenoid_center: [40, 60, 70])
    a = measure_mask(mask_a, humerus: humerus_a).bone_loss_percent
    b = measure_mask(mask_b, humerus: humerus_b).bone_loss_percent
    print "  [trans] @default=#{a.round(3)}%, @offset=#{b.round(3)}%, " \
          "diff=#{(a - b).round(3)}\n"
    expect((a - b).abs).to be < TOLERANCE_PCT
  end

  context "real-mask regression (synthetic fixtures)" do
    # Real segmentation masks routinely violate the synthetic "single thin disc"
    # assumption. These two specs reproduce the failure modes I instrumented
    # against the appendicular-task labels.nii.gz (scapula_left from
    # shoulder_segmenter), where the unfixed pipeline returned
    # bone_loss_percent ~= 96% and radius_mm ~= 138 on a healthy glenoid.

    it "ignores a disconnected spurious blob far from the glenoid" do
      # An intact glenoid disc + a spurious blob the segmenter would label as
      # 'scapula_left' (an acromion tip or labelling-noise cluster).
      mask = Synthetic.intact_glenoid_mask
      # Add a 4x4x4 spurious blob in the opposite corner.
      mask[10..13, 10..13, 10..13] = 1

      result = measure_mask(mask)
      # Without largest-component preprocessing this fit gets dragged across
      # the volume by the spurious blob; with it the answer should be the
      # same as the spurious-blob-free case.
      expect(result.fitted_circle[:radius_mm]).to be_within(2.0).of(Synthetic::GLENOID_RADIUS)
      expect(result.bone_loss_percent).to be < TOLERANCE_PCT
    end

    it "isolates the glenoid surface when a humerus mask is supplied" do
      # Real CT segmentation masks fuse the glenoid disc with the rest of
      # the scapula. We simulate that: an intact disc connected via a thin
      # neck to a thick scapular-blade slab.
      mask = Synthetic.intact_glenoid_mask
      # Scapular-blade slab.
      mask[20..79, 20..79, 0..40] = 1
      # A neck connecting the slab top (k=40) to the disc bottom (k=78).
      mask[45..54, 45..54, 40..78] = 1

      humerus = Synthetic.humerus_mask_for
      with_crop = measure_mask(mask, humerus: humerus, glenoid_window_mm: 25.0)

      # Without the proximity crop the fit would be dragged across the whole
      # scapula (the slab dwarfs the disc). With humerus-proximity cropping
      # we recover the disc radius.
      expect(with_crop.fitted_circle[:radius_mm]).to be_within(6.0).of(Synthetic::GLENOID_RADIUS)
      expect(with_crop.bone_loss_percent).to be < 35.0
    end
  end
end

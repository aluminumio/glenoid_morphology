# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::SurfaceExtraction do
  it "finds the surface of a solid cube" do
    mask = Numo::Bit.zeros(10, 10, 10)
    mask[2..7, 2..7, 2..7] = 1
    # Interior is 4x4x4 = 64 voxels, total cube 6x6x6 = 216, so surface = 152.
    voxels = described_class.surface_voxels(mask)
    expect(voxels.shape[0]).to eq(216 - 64)
  end

  it "returns empty array for empty mask" do
    mask = Numo::Bit.zeros(8, 8, 8)
    voxels = described_class.surface_voxels(mask)
    expect(voxels.shape[0]).to eq(0)
  end

  it "every interior voxel of a thick blob is excluded" do
    mask = Numo::Bit.zeros(20, 20, 20)
    mask[5..14, 5..14, 5..14] = 1
    surf = described_class.surface_voxels(mask)
    # An interior point [8,8,8] must not appear in surf.
    set = surf.to_a.map { |row| row }
    expect(set).not_to include([8, 8, 8])
    expect(set).to include([5, 5, 5]) # corner is on the surface
  end

  it "filter_facing keeps only points on one side of the centroid" do
    mask = Numo::Bit.zeros(20, 20, 20)
    mask[5..14, 5..14, 5..14] = 1
    surf = described_class.surface_voxels(mask)
    facing = described_class.filter_facing(surf, direction: [0.0, 0.0, 1.0])
    # Centroid k ~ 9.5; facing should keep all points with k > 9.5.
    centroid_k = (Numo::DFloat.cast(surf)[true, 2].mean)
    facing_ks = facing[true, 2].to_a
    expect(facing_ks.all? { |k| k > centroid_k - 0.001 }).to be(true)
    expect(facing.shape[0]).to be < surf.shape[0]
  end

  describe ".largest_component" do
    it "keeps a single blob unchanged" do
      mask = Numo::Bit.zeros(10, 10, 10)
      mask[2..7, 2..7, 2..7] = 1
      cc = described_class.largest_component(mask)
      expect(cc.count_true).to eq(mask.count_true)
    end

    it "drops disconnected spurious blobs" do
      mask = Numo::Bit.zeros(20, 20, 20)
      # Big blob (the "scapula")
      mask[5..14, 5..14, 5..14] = 1
      # Two stray blobs (the "acromion tip" + "noise")
      mask[18..19, 18..19, 18..19] = 1
      mask[0..1, 0..1, 0..1] = 1

      cc = described_class.largest_component(mask)
      expect(cc[10, 10, 10]).to eq(1)
      expect(cc[18, 18, 18]).to eq(0)
      expect(cc[0, 0, 0]).to eq(0)
      expect(cc.count_true).to eq(10 * 10 * 10)
    end

    it "returns an all-zero mask for empty input" do
      mask = Numo::Bit.zeros(5, 5, 5)
      cc = described_class.largest_component(mask)
      expect(cc.count_true).to eq(0)
    end
  end

  describe ".filter_near_anchor" do
    it "keeps voxels within d_min + window of the anchor and drops the rest" do
      idx = Numo::Int64[[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 100]]
      mm  = Numo::DFloat.cast(idx) # 1 mm spacing
      anchor = Numo::DFloat[0.0, 0.0, 0.0]
      kept = described_class.filter_near_anchor(idx, voxels_mm: mm, anchor_mm: anchor, window_mm: 1.5)
      # d_min = 0; window = 1.5 → keep voxels at distance <= 1.5
      # That's [0,0,0] (d=0), [0,0,1] (d=1). Not [0,0,2] (d=2), not [0,0,100].
      expect(kept.shape[0]).to eq(2)
    end
  end
end

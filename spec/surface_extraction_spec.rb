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
end

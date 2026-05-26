# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::Affine do
  it "identity affine returns the input unchanged" do
    a = described_class.identity(spacing_mm: 1.0)
    pts = Numo::Int64[[1, 2, 3], [4, 5, 6]]
    mm = described_class.voxels_to_mm(pts, a)
    expect(mm.to_a).to eq([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  end

  it "isotropic spacing scales coordinates" do
    a = described_class.identity(spacing_mm: 0.5)
    pts = Numo::Int64[[2, 4, 6]]
    mm = described_class.voxels_to_mm(pts, a)
    expect(mm[0, 0]).to be_within(1e-9).of(1.0)
    expect(mm[0, 1]).to be_within(1e-9).of(2.0)
    expect(mm[0, 2]).to be_within(1e-9).of(3.0)
    expect(described_class.voxel_spacing_mm(a)).to be_within(1e-9).of(0.5)
  end

  it "accepts plain ruby array affine" do
    a = [
      [1, 0, 0, 10],
      [0, 1, 0, 20],
      [0, 0, 1, 30],
      [0, 0, 0, 1]
    ]
    mm = described_class.voxels_to_mm(Numo::Int64[[0, 0, 0]], a)
    expect(mm.to_a).to eq([[10.0, 20.0, 30.0]])
  end
end

# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::Projector do
  it "project + lift is an identity for in-plane points" do
    proj = described_class.new(origin: [1.0, 2.0, 3.0], normal: [0.0, 0.0, 1.0])
    pts3d = Numo::DFloat[[2.0, 4.0, 3.0], [-1.0, 0.0, 3.0], [5.0, 5.0, 3.0]]
    pts2d = proj.project(pts3d)
    back  = proj.lift(pts2d)
    err = ((pts3d - back) ** 2).sum
    expect(err).to be_within(1e-10).of(0.0)
  end

  it "drops the normal component for off-plane points" do
    proj = described_class.new(origin: [0.0, 0.0, 0.0], normal: [0.0, 0.0, 1.0])
    pts = Numo::DFloat[[2.0, 4.0, 99.0]]
    proj2 = proj.project(pts)
    expect(proj2[0, 0]).to be_within(1e-9).of(2.0)
    expect(proj2[0, 1]).to be_within(1e-9).of(4.0)
  end

  it "works with a tilted plane" do
    n = Math.sqrt(3)
    proj = described_class.new(origin: [0.0, 0.0, 0.0], normal: [1.0 / n, 1.0 / n, 1.0 / n])
    # A point in the plane: (1, -1, 0) satisfies x+y+z = 0
    pts = Numo::DFloat[[1.0, -1.0, 0.0]]
    pts2d = proj.project(pts)
    back  = proj.lift(pts2d)
    err = ((pts - back) ** 2).sum
    expect(err).to be_within(1e-10).of(0.0)
  end
end

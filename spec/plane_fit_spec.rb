# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::PlaneFit do
  def synthesise_plane_points(normal:, offset:, n: 400, noise_sd: 0.05, seed: 42)
    rng = Random.new(seed)
    nrm = normal.dup
    norm = Math.sqrt(nrm.map { |x| x * x }.sum)
    nrm = nrm.map { |x| x / norm }
    # Build two in-plane basis vectors
    seed_axis = (nrm[0].abs > 0.9 ? [0.0, 1.0, 0.0] : [1.0, 0.0, 0.0])
    dot = nrm.zip(seed_axis).sum { |a, b| a * b }
    u = seed_axis.zip(nrm).map { |s, p| s - dot * p }
    un = Math.sqrt(u.map { |x| x * x }.sum)
    u.map! { |x| x / un }
    v = [
      nrm[1] * u[2] - nrm[2] * u[1],
      nrm[2] * u[0] - nrm[0] * u[2],
      nrm[0] * u[1] - nrm[1] * u[0]
    ]
    p = Numo::DFloat.zeros(n, 3)
    n.times do |i|
      a = (rng.rand - 0.5) * 20
      b = (rng.rand - 0.5) * 20
      noise = (rng.rand - 0.5) * 2 * noise_sd
      pt = [0, 1, 2].map do |k|
        a * u[k] + b * v[k] + (offset + noise) * nrm[k]
      end
      p[i, 0] = pt[0]
      p[i, 1] = pt[1]
      p[i, 2] = pt[2]
    end
    p
  end

  it "recovers the plane normal of a noisy planar point cloud" do
    points = synthesise_plane_points(normal: [0.0, 0.0, 1.0], offset: 3.0)
    plane = described_class.fit(points)
    # Normal should be parallel or anti-parallel to (0,0,1)
    n = plane.normal
    expect(n[2].abs).to be > 0.99
    expect(n[0].abs).to be < 0.05
    expect(n[1].abs).to be < 0.05
    expect(plane.d.abs).to be_within(0.2).of(3.0)
    expect(plane.residual_rms).to be < 0.1
  end

  it "recovers a tilted plane" do
    nrm = [1.0, 1.0, 1.0]
    norm = Math.sqrt(3)
    nrm.map! { |x| x / norm }
    points = synthesise_plane_points(normal: nrm, offset: 2.5)
    plane = described_class.fit(points)
    cos_sim = plane.normal.zip(nrm).sum { |a, b| a * b }.abs
    expect(cos_sim).to be > 0.99
  end

  it "raises on degenerate input" do
    p = Numo::DFloat[[0, 0, 0], [1, 1, 1]]
    expect { described_class.fit(p) }.to raise_error(ArgumentError)
  end
end

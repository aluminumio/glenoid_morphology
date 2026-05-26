# frozen_string_literal: true

require "spec_helper"

RSpec.describe GlenoidMorphology::CircleFit do
  def synth_circle(cx:, cy:, r:, n: 200, noise: 0.05, gap_deg: 0.0, seed: 1)
    rng = Random.new(seed)
    angle_used = (360.0 - gap_deg) * Math::PI / 180.0
    p = Numo::DFloat.zeros(n, 2)
    n.times do |i|
      theta = (i.to_f / (n - 1)) * angle_used
      noise_r = (rng.rand - 0.5) * 2 * noise
      p[i, 0] = cx + (r + noise_r) * Math.cos(theta)
      p[i, 1] = cy + (r + noise_r) * Math.sin(theta)
    end
    p
  end

  it "recovers a centred unit circle (Kasa)" do
    points = synth_circle(cx: 0.0, cy: 0.0, r: 1.0, noise: 0.005)
    c = described_class.kasa(points)
    expect(c.cx).to be_within(0.01).of(0.0)
    expect(c.cy).to be_within(0.01).of(0.0)
    expect(c.radius).to be_within(0.01).of(1.0)
    expect(c.rms_residual).to be < 0.01
  end

  it "recovers an off-centre larger circle (Kasa)" do
    points = synth_circle(cx: 12.5, cy: -3.2, r: 8.0, noise: 0.1, n: 300)
    c = described_class.kasa(points)
    expect(c.cx).to be_within(0.2).of(12.5)
    expect(c.cy).to be_within(0.2).of(-3.2)
    expect(c.radius).to be_within(0.2).of(8.0)
  end

  it "RANSAC is robust to outliers" do
    clean = synth_circle(cx: 5.0, cy: 5.0, r: 10.0, n: 200, noise: 0.1, seed: 7)
    rng = Random.new(13)
    # Add 40 outliers in a totally different region
    outliers = Numo::DFloat.zeros(40, 2)
    40.times do |i|
      outliers[i, 0] = (rng.rand - 0.5) * 4 + 40
      outliers[i, 1] = (rng.rand - 0.5) * 4 + 40
    end
    points = Numo::DFloat.zeros(240, 2)
    points[0...200, true] = clean
    points[200..-1, true] = outliers

    c = described_class.ransac(points, threshold: 1.0, iterations: 400, seed: 99)
    expect(c.cx).to be_within(0.5).of(5.0)
    expect(c.cy).to be_within(0.5).of(5.0)
    expect(c.radius).to be_within(0.5).of(10.0)
  end

  it "fits even when input has a gap (partial arc)" do
    points = synth_circle(cx: 0.0, cy: 0.0, r: 1.0, n: 200, noise: 0.02, gap_deg: 60.0, seed: 5)
    c = described_class.kasa(points)
    expect(c.radius).to be_within(0.05).of(1.0)
  end
end

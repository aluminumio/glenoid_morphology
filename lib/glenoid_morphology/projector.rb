# frozen_string_literal: true

require "numo/narray"

module GlenoidMorphology
  # Project 3D points onto a 2D coordinate system on a plane, and lift 2D
  # points back into 3D.
  #
  # Given a plane (centroid + normal), we pick an orthonormal in-plane basis
  # (u, v) deterministically: u is the projection of the world x-axis onto
  # the plane (falling back to y-axis if x is parallel to the normal),
  # v = normal x u.
  class Projector
    attr_reader :origin, :u, :v, :normal

    # @param plane [PlaneFit::Plane]
    def self.from_plane(plane)
      new(origin: plane.centroid, normal: plane.normal)
    end

    def initialize(origin:, normal:)
      @origin = Numo::DFloat.cast(origin)
      @normal = unit(Numo::DFloat.cast(normal))
      @u, @v = build_basis(@normal)
    end

    # Project [N, 3] 3D points onto the plane, returning [N, 2] (u, v) coords.
    def project(points)
      p = Numo::DFloat.cast(points)
      rel = p - @origin.reshape(1, 3)
      uu = rel.dot(@u)
      vv = rel.dot(@v)
      stack(uu, vv)
    end

    # Lift [N, 2] (u, v) coords back into 3D world coords.
    def lift(points2d)
      p = Numo::DFloat.cast(points2d)
      n = p.shape[0]
      out = Numo::DFloat.zeros(n, 3)
      out[true, 0] = @origin[0] + p[true, 0] * @u[0] + p[true, 1] * @v[0]
      out[true, 1] = @origin[1] + p[true, 0] * @u[1] + p[true, 1] * @v[1]
      out[true, 2] = @origin[2] + p[true, 0] * @u[2] + p[true, 1] * @v[2]
      out
    end

    private

    def unit(v)
      v / Math.sqrt((v ** 2).sum)
    end

    def build_basis(n)
      x_axis = Numo::DFloat[1.0, 0.0, 0.0]
      seed = (n.dot(x_axis).abs > 0.9) ? Numo::DFloat[0.0, 1.0, 0.0] : x_axis
      u = seed - n * n.dot(seed)
      u = unit(u)
      v = cross(n, u)
      v = unit(v)
      [u, v]
    end

    def cross(a, b)
      Numo::DFloat[
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
      ]
    end

    def stack(uu, vv)
      n = uu.size
      out = Numo::DFloat.zeros(n, 2)
      out[true, 0] = uu
      out[true, 1] = vv
      out
    end
  end
end

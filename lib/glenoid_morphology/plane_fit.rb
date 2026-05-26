# frozen_string_literal: true

require "numo/narray"
require "numo/linalg"

module GlenoidMorphology
  # PCA-based plane fit to a 3D point cloud.
  #
  # The plane normal is the eigenvector of the covariance matrix with the
  # smallest eigenvalue; the plane passes through the centroid.
  #
  # Returned as (a, b, c, d) such that a*x + b*y + c*z = d, with (a, b, c)
  # a unit vector and d the signed offset from origin.
  module PlaneFit
    Plane = Struct.new(:a, :b, :c, :d, :centroid, :normal, :residual_rms) do
      def coefficients
        [a, b, c, d]
      end

      # Signed distance of each [N,3] point to the plane.
      def signed_distance(points)
        p = Numo::DFloat.cast(points)
        n = Numo::DFloat.cast([a, b, c]).reshape(3, 1)
        (p.dot(n)).reshape(p.shape[0]) - d
      end
    end

    module_function

    # @param points [Numo::DFloat] [N, 3]
    # @return [Plane]
    def fit(points)
      p = Numo::DFloat.cast(points)
      raise ArgumentError, "need at least 3 points" if p.shape[0] < 3

      centroid = p.mean(axis: 0)
      centred  = p - centroid.reshape(1, 3)
      cov      = centred.transpose.dot(centred) / p.shape[0].to_f

      # eigh returns ascending eigenvalues; smallest first.
      eigvals, eigvecs = Numo::Linalg.eigh(cov)
      normal = eigvecs[true, 0].dup
      # Normalise (defensive — eigh returns unit eigenvectors but float drift is real).
      normal /= Math.sqrt((normal ** 2).sum)

      d = normal.dot(centroid)
      residual_rms = Math.sqrt(eigvals[0].abs) # smallest eigenvalue == variance along normal

      Plane.new(
        normal[0], normal[1], normal[2], d,
        [centroid[0], centroid[1], centroid[2]],
        [normal[0], normal[1], normal[2]],
        residual_rms
      )
    end
  end
end

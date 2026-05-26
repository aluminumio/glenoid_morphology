# frozen_string_literal: true

require "numo/narray"
require "numo/linalg"

module GlenoidMorphology
  # 2D circle fitting via Kasa's algebraic least-squares method,
  # optionally wrapped in RANSAC to reject outliers.
  #
  # Kasa minimises sum_i (x_i^2 + y_i^2 - 2*a*x_i - 2*b*y_i - c)^2, where
  # the circle is (x-a)^2 + (y-b)^2 = r^2 with r^2 = a^2 + b^2 + c.
  module CircleFit
    Circle = Struct.new(:cx, :cy, :radius, :rms_residual, :inlier_ratio) do
      def center
        [cx, cy]
      end

      def signed_radial_residuals(points)
        p = Numo::DFloat.cast(points)
        dx = p[true, 0] - cx
        dy = p[true, 1] - cy
        Numo::NMath.sqrt(dx * dx + dy * dy) - radius
      end
    end

    module_function

    # Direct Kasa-method least-squares fit.
    # @param points [Numo::DFloat] [N, 2]
    # @return [Circle]
    def kasa(points)
      p = Numo::DFloat.cast(points)
      n = p.shape[0]
      raise ArgumentError, "need at least 3 points" if n < 3

      x = p[true, 0]
      y = p[true, 1]

      a_mat = Numo::DFloat.zeros(n, 3)
      a_mat[true, 0] = 2.0 * x
      a_mat[true, 1] = 2.0 * y
      a_mat[true, 2] = 1.0
      b_vec = x * x + y * y

      sol, = Numo::Linalg.lstsq(a_mat, b_vec)
      cx = sol[0]
      cy = sol[1]
      c  = sol[2]
      r2 = c + cx * cx + cy * cy
      r  = Math.sqrt([r2, 0.0].max)

      residuals = (Numo::NMath.sqrt((x - cx) ** 2 + (y - cy) ** 2) - r)
      rms = Math.sqrt((residuals ** 2).mean)
      Circle.new(cx, cy, r, rms, 1.0)
    end

    # RANSAC + Kasa refinement.
    #
    # @param points [Numo::DFloat] [N, 2]
    # @param threshold [Float] inlier radial-distance tolerance (mm)
    # @param iterations [Integer]
    # @param min_inlier_ratio [Float]
    # @param seed [Integer, nil] RNG seed for reproducible specs
    # @return [Circle]
    def ransac(points, threshold: 1.5, iterations: 200, min_inlier_ratio: 0.4, seed: nil)
      p = Numo::DFloat.cast(points)
      n = p.shape[0]
      raise ArgumentError, "need at least 3 points" if n < 3

      rng = seed ? Random.new(seed) : Random.new
      best_inliers = nil
      best_score = -1

      iterations.times do
        idx = sample3(n, rng)
        sample = p[idx, true]
        next unless sample.shape[0] == 3

        circle = three_point_circle(sample)
        next unless circle
        next if circle.radius <= 0 || !circle.radius.finite?

        residuals = (Numo::NMath.sqrt((p[true, 0] - circle.cx) ** 2 +
                                      (p[true, 1] - circle.cy) ** 2) - circle.radius).abs
        inlier_mask = residuals < threshold
        count = inlier_mask.count_true

        if count > best_score
          best_score = count
          best_inliers = inlier_mask
        end
      end

      if best_inliers.nil? || best_score.to_f / n < min_inlier_ratio
        return kasa(p) # graceful fallback
      end

      inlier_idx = best_inliers.where
      refined = kasa(p[inlier_idx, true])
      Circle.new(refined.cx, refined.cy, refined.radius, refined.rms_residual,
                 best_score.to_f / n)
    end

    # ---- helpers ----

    def sample3(n, rng)
      raise ArgumentError, "need at least 3 points" if n < 3

      a = rng.rand(n)
      b = a
      b = rng.rand(n) while b == a
      c = a
      c = rng.rand(n) while c == a || c == b
      Numo::Int32[a, b, c]
    end

    # Exact circle through three points. Returns nil if collinear.
    def three_point_circle(points)
      x1 = points[0, 0]; y1 = points[0, 1]
      x2 = points[1, 0]; y2 = points[1, 1]
      x3 = points[2, 0]; y3 = points[2, 1]

      ax = x2 - x1; ay = y2 - y1
      bx = x3 - x1; by = y3 - y1
      d = 2.0 * (ax * by - ay * bx)
      return nil if d.abs < 1e-12

      a_sq = ax * ax + ay * ay
      b_sq = bx * bx + by * by
      ux = (by * a_sq - ay * b_sq) / d
      uy = (ax * b_sq - bx * a_sq) / d

      cx = x1 + ux
      cy = y1 + uy
      r  = Math.sqrt(ux * ux + uy * uy)
      Circle.new(cx, cy, r, 0.0, 0.0)
    end
  end
end

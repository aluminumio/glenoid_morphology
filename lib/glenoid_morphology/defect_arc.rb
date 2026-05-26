# frozen_string_literal: true

require "numo/narray"

module GlenoidMorphology
  # Given a best-fit circle in 2D and the set of 2D rim points, find the
  # angular arc(s) that lack supporting rim points — that is the bone defect.
  #
  # We bin the 2D rim points into angular slots around the circle centre and
  # mark each slot as "bone present" if at least one point is within
  # `radial_tolerance` of the circle radius. Contiguous empty slots are the
  # candidate defects; we return the LARGEST one as the principal defect.
  module DefectArc
    Arc = Struct.new(:start_angle, :end_angle, :arc_length_mm, :bone_loss_fraction,
                     :area_lost_mm2, :area_total_mm2, :angular_span_rad,
                     :occupancy) do
      def bone_loss_percent
        bone_loss_fraction * 100.0
      end
    end

    module_function

    # @param points2d [Numo::DFloat] [N, 2] points on the rim plane
    # @param circle [CircleFit::Circle] fitted circle
    # @param radial_tolerance [Float] mm. How far from the circle edge counts as "on the rim"
    # @param bins [Integer] number of angular bins (default 360 == 1 deg)
    # @return [Arc]
    def find(points2d, circle:, radial_tolerance: 2.0, bins: 360)
      p = Numo::DFloat.cast(points2d)
      n = p.shape[0]
      occupancy = Numo::Bit.zeros(bins)

      if n.positive?
        dx = p[true, 0] - circle.cx
        dy = p[true, 1] - circle.cy
        radii = Numo::NMath.sqrt(dx * dx + dy * dy)
        angles = Numo::DFloat.zeros(n)
        n.times { |i| angles[i] = Math.atan2(dy[i], dx[i]) } # atan2 not vectorised in Numo
        # Normalise to [0, 2pi)
        twopi = 2.0 * Math::PI
        angles = ((angles + twopi) % twopi)

        radial_residual = (radii - circle.radius).abs
        on_rim = radial_residual < radial_tolerance

        on_rim.where.to_a.each do |idx|
          bin = (angles[idx] / twopi * bins).floor % bins
          occupancy[bin] = 1
        end
      end

      gap = largest_gap(occupancy, bins)
      angular_span = gap[:span] * (2.0 * Math::PI / bins)
      arc_length   = angular_span * circle.radius

      # Circular segment area (the "bone lost" area):
      # A_segment = 0.5 * r^2 * (theta - sin(theta))
      area_lost = 0.5 * (circle.radius**2) * (angular_span - Math.sin(angular_span))
      area_total = Math::PI * (circle.radius**2)
      fraction = area_total.positive? ? (area_lost / area_total) : 0.0

      start_angle = gap[:start] * (2.0 * Math::PI / bins)
      end_angle   = gap[:end_excl] * (2.0 * Math::PI / bins)

      Arc.new(start_angle, end_angle, arc_length, fraction,
              area_lost, area_total, angular_span, occupancy)
    end

    # Find the largest run of zeros in a circular occupancy bitmap.
    # Returns { start:, end_excl:, span: } in bin units.
    # If the whole map is zero we return the full circle.
    # If the whole map is one we return zero span.
    def largest_gap(occupancy, bins)
      occ = occupancy.to_a
      if occ.all? { |b| b.zero? }
        return { start: 0, end_excl: bins, span: bins }
      end
      if occ.all? { |b| b == 1 }
        return { start: 0, end_excl: 0, span: 0 }
      end

      # Double the array to handle wrap-around, scan for longest run of zeros
      # whose length is at most `bins` (anything longer is degenerate).
      doubled = occ + occ
      best_len = 0
      best_start = 0
      i = 0
      while i < doubled.length
        if doubled[i].zero?
          j = i
          j += 1 while j < doubled.length && doubled[j].zero?
          run = j - i
          if run > best_len
            best_len = [run, bins].min
            best_start = i % bins
          end
          i = j
        else
          i += 1
        end
      end

      { start: best_start, end_excl: (best_start + best_len) % bins, span: best_len }
    end
  end
end

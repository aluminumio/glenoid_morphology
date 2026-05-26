# frozen_string_literal: true

require "numo/linalg"
require_relative "affine"
require_relative "surface_extraction"
require_relative "plane_fit"
require_relative "projector"
require_relative "circle_fit"
require_relative "defect_arc"
require_relative "measurement"

module GlenoidMorphology
  # Top-level orchestration. See GlenoidMorphology.measure.
  module Measure
    DEFAULT_OPTS = {
      humerus_mask: nil,
      side: nil,                # :left, :right, or nil to infer
      ransac_iterations: 200,
      ransac_threshold_mm: 1.5,
      rim_radial_tolerance_mm: 1.0,
      defect_bins: 360,
      seed: nil
    }.freeze

    module_function

    def call(scapula_mask:, affine:, **opts)
      opts = DEFAULT_OPTS.merge(opts)

      # 1. Surface voxels of the scapula (in voxel index space).
      surface_idx = SurfaceExtraction.surface_voxels(scapula_mask)
      raise Error, "scapula mask has no surface voxels" if surface_idx.shape[0].zero?

      # 2. Filter to the glenoid-facing side.
      facing_dir = infer_glenoid_direction(scapula_mask,
                                           surface_idx: surface_idx,
                                           humerus_mask: opts[:humerus_mask],
                                           side: opts[:side])
      glenoid_idx = SurfaceExtraction.filter_facing(surface_idx, direction: facing_dir)
      glenoid_idx = surface_idx if glenoid_idx.shape[0] < 10

      # 3. Convert voxel indices to mm.
      surface_mm = Affine.voxels_to_mm(glenoid_idx, affine)

      # 4. En-face plane via PCA.
      plane = PlaneFit.fit(surface_mm)

      # 5. Project surface points onto the plane.
      projector = Projector.from_plane(plane)
      points2d  = projector.project(surface_mm)

      # 6. Per-angle outermost rim points (one point per angular bin).
      rim_points2d = select_rim_candidates(points2d, bins: opts[:defect_bins])

      # 7. Fit a circle. The rim radii distribute bimodally: arc points at the
      #    true rim radius and chord-defect points at a smaller radius. We
      #    seed with the top tertile, but then run a RANSAC pass on the full
      #    rim — RANSAC's 3-point sample inherently finds the dominant circle
      #    (the arc) since the arc has more points than the chord — and finally
      #    iteratively reweight to converge on the largest-inlier-set circle.
      fit_seed = top_quantile_by_radius(rim_points2d, quantile: 0.35)
      seed_circle = CircleFit.ransac(
        fit_seed,
        threshold: opts[:ransac_threshold_mm],
        iterations: opts[:ransac_iterations],
        seed: opts[:seed]
      )
      circle = CircleFit.ransac(
        rim_points2d,
        threshold: opts[:ransac_threshold_mm],
        iterations: opts[:ransac_iterations] * 3,
        seed: opts[:seed]
      )
      # Pick whichever has the larger radius (the arc-only fit is biased small
      # by chord-cut points that crept in; whichever fit gives the larger
      # radius is closer to the intact rim).
      circle = seed_circle if seed_circle.radius > circle.radius
      circle = refine_circle(rim_points2d, circle, tolerance: opts[:rim_radial_tolerance_mm])

      # 8. Defect arc. Walk the fitted circle in 2D; for each angular bin,
      #    sample the *mask itself* (in voxel space) along a short ray inward
      #    from the rim point to ask "is there bone here?". An angle without
      #    bone within tolerance is a defect bin.
      arc = walk_circle_for_defect(
        circle: circle,
        projector: projector,
        affine: affine,
        scapula_mask: scapula_mask,
        bins: opts[:defect_bins],
        tolerance_mm: opts[:rim_radial_tolerance_mm]
      )

      # 9. Lift rim & centre back to 3D mm for output.
      rim_3d  = projector.lift(rim_points2d).to_a
      center3d = projector.lift(Numo::DFloat[[circle.cx, circle.cy]])[0, true].to_a

      Measurement.new(
        en_face_plane: { a: plane.a, b: plane.b, c: plane.c, d: plane.d },
        glenoid_rim_points: rim_3d,
        fitted_circle: {
          center: center3d,
          radius_mm: circle.radius,
          normal: plane.normal
        },
        defect_arc: {
          start_angle: arc.start_angle,
          end_angle: arc.end_angle,
          arc_length_mm: arc.arc_length_mm
        },
        bone_loss_percent: arc.bone_loss_percent,
        confidence: confidence_from(circle: circle, plane: plane, n_points: rim_points2d.shape[0])
      )
    end

    # Outward direction (in voxel index space) from the scapula centroid
    # toward where the glenoid should be. Strategy:
    #   1. If a humerus mask is given, take centroid_humerus - centroid_scapula.
    #   2. Otherwise infer from `side`: anatomically the glenoid faces laterally.
    #      With no body-axis info we approximate by picking the axis along which
    #      the scapula has the largest extent, then point AWAY from the bulk.
    def infer_glenoid_direction(scapula_mask, surface_idx:, humerus_mask:, side:)
      if humerus_mask
        h_idx = SurfaceExtraction.indices_of_true(humerus_mask.cast_to(Numo::Bit))
        if h_idx.shape[0].positive?
          h_centroid = Numo::DFloat.cast(h_idx).mean(axis: 0)
          s_centroid = Numo::DFloat.cast(surface_idx).mean(axis: 0)
          return (h_centroid - s_centroid).to_a
        end
      end

      # Fallback: pick the axis with the smallest scapula extent (the body is
      # roughly planar — the smallest-extent axis is the through-the-blade
      # direction, which is also approximately the glenoid-facing axis).
      v = Numo::DFloat.cast(surface_idx)
      span = v.max(axis: 0) - v.min(axis: 0)
      axis = span.to_a.each_with_index.min_by { |a, _| a }[1]
      dir = [0.0, 0.0, 0.0]
      dir[axis] = (side == :left ? -1.0 : 1.0)
      dir
    end

    # Extract rim points: the convex hull of the projected 2D point cloud.
    # The articular surface of a glenoid (intact OR chord-defect) is convex,
    # so its rim is exactly the convex hull boundary of the projection.
    # This sidesteps the bimodal-radius problem of a per-angle max-r approach.
    def select_rim_candidates(points2d, bins: nil) # rubocop:disable Lint/UnusedMethodArgument
      n = points2d.shape[0]
      return points2d if n < 4

      pts = points2d.to_a
      hull = convex_hull(pts)
      out = Numo::DFloat.zeros(hull.size, 2)
      hull.each_with_index do |(x, y), i|
        out[i, 0] = x
        out[i, 1] = y
      end
      out
    end

    # Andrew's monotone-chain convex hull in O(n log n). Returns CCW vertices.
    def convex_hull(points)
      pts = points.uniq.sort
      return pts if pts.size <= 2

      cross = ->(o, a, b) { (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]) }
      lower = []
      pts.each do |p|
        lower.pop while lower.size >= 2 && cross.call(lower[-2], lower[-1], p) <= 0
        lower << p
      end
      upper = []
      pts.reverse_each do |p|
        upper.pop while upper.size >= 2 && cross.call(upper[-2], upper[-1], p) <= 0
        upper << p
      end
      lower[0...-1] + upper[0...-1]
    end

    # Walk the fitted circle in the rim-plane 2D frame. For each angular bin,
    # cast a ray FROM the centre OUTWARD past the fitted rim, and record the
    # maximum radius at which bone is still present. A bin is "intact" iff
    # max-bone-radius >= fitted_radius - tolerance; otherwise it's a defect bin
    # (bone is recessed inward of the rim, which is exactly how glenoid bone
    # loss manifests).
    def walk_circle_for_defect(circle:, projector:, affine:, scapula_mask:,
                               bins:, tolerance_mm:)
      occupancy = Numo::Bit.zeros(bins)
      affine_n = Affine.to_narray(affine)
      inv = Numo::Linalg.inv(affine_n)
      shape = scapula_mask.shape

      # Sample radii from 0 out to radius + tolerance, step ~ half a voxel.
      voxel_mm = Affine.voxel_spacing_mm(affine_n)
      step = [voxel_mm * 0.5, 0.25].max
      r_max = circle.radius + tolerance_mm
      n_samples = (r_max / step).ceil + 1

      # Through-plane offsets: a few samples ALONG the plane normal centred on
      # the in-plane sample point. Compensates for tilt between the fitted
      # plane and the true disc plane.
      normal = projector.normal
      n_through = 5
      through_spacing = voxel_mm

      bins.times do |b|
        theta = (b.to_f / bins) * 2.0 * Math::PI
        cos_t = Math.cos(theta)
        sin_t = Math.sin(theta)

        max_bone_r = -1.0
        n_samples.times do |s|
          r_sample = s * step
          u = circle.cx + r_sample * cos_t
          v = circle.cy + r_sample * sin_t
          point3d = projector.lift(Numo::DFloat[[u, v]])[0, true]

          hit = false
          (-(n_through / 2)..(n_through / 2)).each do |t|
            off = t * through_spacing
            p3 = Numo::DFloat[
              point3d[0] + off * normal[0],
              point3d[1] + off * normal[1],
              point3d[2] + off * normal[2],
              1.0
            ]
            vox = inv.dot(p3)
            i = vox[0].round; j = vox[1].round; k = vox[2].round
            next if i < 0 || j < 0 || k < 0 || i >= shape[0] || j >= shape[1] || k >= shape[2]
            if scapula_mask[i, j, k] == 1
              hit = true
              break
            end
          end
          max_bone_r = r_sample if hit
        end
        occupancy[b] = 1 if max_bone_r >= circle.radius - tolerance_mm
      end

      gap = DefectArc.largest_gap(occupancy, bins)
      angular_span = gap[:span] * (2.0 * Math::PI / bins)
      arc_length = angular_span * circle.radius
      area_lost = 0.5 * (circle.radius**2) * (angular_span - Math.sin(angular_span))
      area_total = Math::PI * (circle.radius**2)
      fraction = area_total.positive? ? (area_lost / area_total) : 0.0

      start_angle = gap[:start] * (2.0 * Math::PI / bins)
      end_angle   = gap[:end_excl] * (2.0 * Math::PI / bins)
      DefectArc::Arc.new(start_angle, end_angle, arc_length, fraction,
                         area_lost, area_total, angular_span, occupancy)
    end

    # Iteratively re-fit the circle using only rim points whose signed radial
    # residual is >= -tolerance (i.e. points ON or OUTSIDE the current circle).
    # Chord-defect points have signed residual << -tolerance and get excluded.
    def refine_circle(rim_points2d, circle, tolerance:, max_iter: 6)
      current = circle
      max_iter.times do
        residuals = current.signed_radial_residuals(rim_points2d)
        keep_mask = residuals >= -tolerance
        kept_idx  = keep_mask.where
        break if kept_idx.size < 8

        kept = rim_points2d[kept_idx, true]
        next_circle = CircleFit.kasa(kept)

        dx = (next_circle.cx - current.cx).abs
        dy = (next_circle.cy - current.cy).abs
        dr = (next_circle.radius - current.radius).abs
        current = next_circle
        break if dx < 0.05 && dy < 0.05 && dr < 0.05
      end
      current
    end

    # Keep only points in the top `quantile` of radial distance from centroid.
    # Used to bias circle-fit toward the intact rim rather than chord-cut points.
    def top_quantile_by_radius(points2d, quantile: 0.5)
      n = points2d.shape[0]
      return points2d if n < 20

      centroid = points2d.mean(axis: 0)
      dx = points2d[true, 0] - centroid[0]
      dy = points2d[true, 1] - centroid[1]
      r = Numo::NMath.sqrt(dx * dx + dy * dy)
      sorted = r.sort
      cutoff = sorted[(n * (1.0 - quantile)).to_i]
      mask = r >= cutoff
      idx = mask.where
      points2d[idx, true]
    end

    def confidence_from(circle:, plane:, n_points:)
      # Heuristic: low circle residual + reasonable inlier ratio + many points => high confidence.
      res_term = 1.0 / (1.0 + circle.rms_residual / [circle.radius * 0.05, 0.1].max)
      inlier_term = circle.inlier_ratio
      pts_term = [n_points / 500.0, 1.0].min
      plane_term = 1.0 / (1.0 + plane.residual_rms / [circle.radius * 0.05, 0.1].max)
      [(res_term * 0.4 + inlier_term * 0.3 + pts_term * 0.1 + plane_term * 0.2), 1.0].min.clamp(0.0, 1.0)
    end
  end
end

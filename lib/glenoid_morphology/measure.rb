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
      method: :pico,            # :pico (Baudi 2005 inferior-2/3 circle) or :full_circle
      ransac_iterations: 200,
      ransac_threshold_mm: 1.5,
      rim_radial_tolerance_mm: 1.0,
      defect_bins: 360,
      # Largest connected component preprocessing. Real segmentations have
      # spurious blobs; disabling for fully synthetic single-blob masks.
      largest_component: true,
      # When a humerus mask is supplied, restrict scapula surface voxels to
      # those within `d_min + glenoid_window_mm` of the humerus centroid.
      # ~30 mm cleanly isolates the glenoid surface from the rest of the
      # scapula at typical CT resolution (real glenoid diameter is ~25 mm).
      glenoid_window_mm: 30.0,
      seed: nil
    }.freeze

    # Pico method constant: keep rim points within +/- 120° (2/3 * pi) of the
    # inferior pole, in the en-face plane's polar coordinate system anchored
    # at the rim centroid with theta=0 at the inferior pole.
    PICO_HALF_ANGLE_RAD = (2.0 / 3.0) * Math::PI

    module_function

    def call(scapula_mask:, affine:, **opts)
      opts = DEFAULT_OPTS.merge(opts)
      method = opts[:method]
      unless %i[pico full_circle].include?(method)
        raise ArgumentError, "method must be :pico or :full_circle, got #{method.inspect}"
      end

      # 0. Strip spurious blobs: keep only the largest connected component.
      #    Skipped if disabled (synthetic single-blob test masks).
      working_mask = opts[:largest_component] ? SurfaceExtraction.largest_component(scapula_mask) : scapula_mask

      # 1. Surface voxels of the scapula (in voxel index space).
      surface_idx = SurfaceExtraction.surface_voxels(working_mask)
      raise Error, "scapula mask has no surface voxels" if surface_idx.shape[0].zero?

      # 2a. If we have a humerus mask, crop surface voxels to the region near
      #     the humerus centroid — this is the glenoid surface. Without this
      #     step the whole scapular blade dominates the PCA plane fit.
      if opts[:humerus_mask]
        hum_idx = SurfaceExtraction.indices_of_true(opts[:humerus_mask].cast_to(Numo::Bit))
        if hum_idx.shape[0].positive?
          hum_centroid_vox = Numo::DFloat.cast(hum_idx).mean(axis: 0)
          hum_centroid_mm  = Affine.voxels_to_mm(hum_centroid_vox.reshape(1, 3), affine)[0, true]
          surface_mm_all   = Affine.voxels_to_mm(surface_idx, affine)
          cropped = SurfaceExtraction.filter_near_anchor(
            surface_idx,
            voxels_mm: surface_mm_all,
            anchor_mm: hum_centroid_mm,
            window_mm: opts[:glenoid_window_mm]
          )
          surface_idx = cropped if cropped.shape[0] >= 30
        end
      end

      # 2b. Filter to the glenoid-facing side (residual cleanup; mostly
      #     redundant when the proximity crop above ran).
      facing_dir = infer_glenoid_direction(working_mask,
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

      # 7. Full-rim circle fit. The rim radii distribute bimodally: arc points at the
      #    true rim radius and chord-defect points at a smaller radius. We
      #    seed with the top tertile, but then run a RANSAC pass on the full
      #    rim — RANSAC's 3-point sample inherently finds the dominant circle
      #    (the arc) since the arc has more points than the chord — and finally
      #    iteratively reweight to converge on the largest-inlier-set circle.
      full_circle = fit_rim_circle(rim_points2d, opts)

      # 8. Optional Pico inferior-2/3 reference circle. Fit only to rim points
      #    in the inferior 240° arc; this segment of a healthy glenoid IS
      #    circular (the pear-shaped narrowing is on the superior 1/3 only).
      pico_circle = nil
      pico_inferior_mask = nil
      si_dir2d = nil
      rim_centroid2d = rim_points2d.mean(axis: 0)
      if method == :pico
        si_dir2d, used_fallback = infer_inferior_direction_2d(
          plane: plane, affine: affine, rim_points2d: rim_points2d
        )
        pico_inferior_mask = inferior_two_thirds_mask(
          rim_points2d, centroid: rim_centroid2d, si_dir2d: si_dir2d
        )
        n_inf = pico_inferior_mask.count_true
        n_total = rim_points2d.shape[0]
        ratio = n_inf.to_f / [n_total, 1].max
        if n_inf < 8 || ratio < 0.20
          raise Error, "Pico inferior-2/3 filter retained too few rim points " \
                       "(#{n_inf}/#{n_total} = #{(ratio * 100).round(1)}%) — " \
                       "SI direction is likely wrong. " \
                       "Try method: :full_circle, or pass a humerus_mask."
        end
        pico_inferior_pts = rim_points2d[pico_inferior_mask.where, true]
        pico_circle = fit_rim_circle(pico_inferior_pts, opts)
        if used_fallback
          warn "[glenoid_morphology] Pico: SI direction is nearly parallel to " \
               "the en-face plane normal; fell back to longest-PCA axis of the " \
               "rim cloud."
        end
      end

      active_circle = method == :pico ? pico_circle : full_circle

      # 9. Defect arc + bone-loss. For Pico, compute fraction over only the
      #    inferior-2/3 sector (the natural pear-narrowing on the superior
      #    1/3 is NOT bone loss and must not be counted).
      pico_sector = method == :pico ? { si_dir2d: si_dir2d, half_angle: PICO_HALF_ANGLE_RAD } : nil
      arc = walk_circle_for_defect(
        circle: active_circle,
        projector: projector,
        affine: affine,
        scapula_mask: working_mask,
        bins: opts[:defect_bins],
        tolerance_mm: opts[:rim_radial_tolerance_mm],
        pico_sector: pico_sector
      )

      # Legacy full-circle bone-loss (always computed, for back-comparison).
      legacy_arc = if method == :full_circle
                     arc
                   else
                     walk_circle_for_defect(
                       circle: full_circle,
                       projector: projector,
                       affine: affine,
                       scapula_mask: working_mask,
                       bins: opts[:defect_bins],
                       tolerance_mm: opts[:rim_radial_tolerance_mm],
                       pico_sector: nil
                     )
                   end

      # 10. Lift rim & circle centres back to 3D mm for output.
      rim_3d   = projector.lift(rim_points2d).to_a
      active_center3d = projector.lift(Numo::DFloat[[active_circle.cx, active_circle.cy]])[0, true].to_a
      pico_hash = if pico_circle
                    pico_center3d = projector.lift(Numo::DFloat[[pico_circle.cx, pico_circle.cy]])[0, true].to_a
                    { center: pico_center3d, radius_mm: pico_circle.radius, normal: plane.normal }
                  end

      Measurement.new(
        en_face_plane: { a: plane.a, b: plane.b, c: plane.c, d: plane.d },
        glenoid_rim_points: rim_3d,
        fitted_circle: {
          center: active_center3d,
          radius_mm: active_circle.radius,
          normal: plane.normal
        },
        pico_circle: pico_hash,
        defect_arc: {
          start_angle: arc.start_angle,
          end_angle: arc.end_angle,
          arc_length_mm: arc.arc_length_mm
        },
        bone_loss_percent: arc.bone_loss_percent,
        legacy_bone_loss_percent: legacy_arc.bone_loss_percent,
        method: method,
        confidence: confidence_from(circle: active_circle, plane: plane, n_points: rim_points2d.shape[0])
      )
    end

    # Wraps the seed-then-refine circle fit used for both the full-rim circle
    # and the Pico inferior-2/3 subset.
    def fit_rim_circle(points2d, opts)
      n = points2d.shape[0]
      raise ArgumentError, "need >= 3 rim points to fit, got #{n}" if n < 3

      seed_pts = top_quantile_by_radius(points2d, quantile: 0.35)
      seed_circle = CircleFit.ransac(
        seed_pts,
        threshold: opts[:ransac_threshold_mm],
        iterations: opts[:ransac_iterations],
        seed: opts[:seed]
      )
      circle = CircleFit.ransac(
        points2d,
        threshold: opts[:ransac_threshold_mm],
        iterations: opts[:ransac_iterations] * 3,
        seed: opts[:seed]
      )
      # Larger-radius fit is closer to the intact arc (chord-cut points bias small).
      circle = seed_circle if seed_circle.radius > circle.radius
      refine_circle(points2d, circle, tolerance: opts[:rim_radial_tolerance_mm])
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

    # Compute the inferior direction expressed in the en-face plane's 2D basis.
    #
    # World-space inferior is the negative of the affine's superior axis. For
    # a standard NIfTI RAS+ affine, "superior" maps to the world +Z direction;
    # the inferior world unit vector is then [0, 0, -1] in mm. (For non-RAS
    # affines, we'd want to pull the column of the affine corresponding to the
    # SI image axis -- but in practice the inputs we consume are RAS.)
    #
    # We project that 3D unit vector onto the en-face plane (subtracting the
    # component along the plane normal), normalise, and dot with (u, v) to get
    # the 2D vector that "points inferior" in the projected coordinate system.
    #
    # @return [Array(Numo::DFloat[2], Boolean)] [si_dir2d, used_fallback]
    def infer_inferior_direction_2d(plane:, affine:, rim_points2d:) # rubocop:disable Lint/UnusedMethodArgument
      world_inf = world_inferior_unit(affine)
      n  = Numo::DFloat.cast(plane.normal)
      n /= Math.sqrt((n ** 2).sum)
      tangential = world_inf - n * world_inf.dot(n)
      norm_t = Math.sqrt((tangential ** 2).sum)

      if norm_t < 0.15
        # Plane normal is nearly parallel to world inferior axis (glenoid would
        # be facing straight up/down — unusual). Fall back to the longest PCA
        # axis of the rim point cloud as the SI proxy.
        return [longest_pca_axis_2d(rim_points2d), true]
      end

      tangential /= norm_t
      # Express in 2D plane basis. Projector picks (u, v) deterministically.
      proj = Projector.from_plane(plane)
      u2 = tangential.dot(proj.u)
      v2 = tangential.dot(proj.v)
      mag = Math.sqrt(u2 * u2 + v2 * v2)
      if mag < 1e-9
        return [longest_pca_axis_2d(rim_points2d), true]
      end

      [Numo::DFloat[u2 / mag, v2 / mag], false]
    end

    # Negative of the affine's superior axis -- i.e. the inferior unit vector
    # in world mm coordinates. For a RAS+ affine the SI image axis is the
    # column most aligned with world Z; inferior is then approximately
    # [0, 0, -1]. We compute it generically by picking the affine column that
    # is most aligned with world Z and flipping if needed -- works for RAS,
    # LPS, and reorientations -- but defaults to [0, 0, -1] if the affine is
    # unavailable.
    def world_inferior_unit(affine)
      a = Affine.to_narray(affine)
      # Each of the first three columns is one image-axis direction in mm.
      best_axis = nil
      best_align = -1.0
      sign = 1.0
      (0..2).each do |i|
        col = Numo::DFloat[a[0, i], a[1, i], a[2, i]]
        nrm = Math.sqrt((col ** 2).sum)
        next if nrm < 1e-9

        unit = col / nrm
        align = unit[2].abs # how aligned with world Z
        if align > best_align
          best_align = align
          best_axis = unit
          sign = unit[2].positive? ? 1.0 : -1.0
        end
      end
      # Inferior is the negative of the superior direction.
      # If best_axis[2] > 0 (axis points superior), inferior = -best_axis.
      best_axis ? best_axis * (-sign) : Numo::DFloat[0.0, 0.0, -1.0]
    end

    # Longest principal axis of a 2D point cloud, returned as a unit 2-vector.
    # Used as the SI fallback when the en-face plane is nearly perpendicular
    # to the world inferior axis.
    def longest_pca_axis_2d(points2d)
      c = points2d.mean(axis: 0)
      centred = points2d - c.reshape(1, 2)
      cov = centred.transpose.dot(centred) / [points2d.shape[0], 1].max.to_f
      _vals, vecs = Numo::Linalg.eigh(cov)
      # eigh returns ascending; the largest-variance axis is the last column.
      axis = vecs[true, -1].dup
      axis /= Math.sqrt((axis ** 2).sum)
      # Sign-disambiguate: point toward the bin with more points (so we have a
      # consistent "inferior pole" direction; either sign technically works
      # for the |theta| <= 120° filter, but consistency is friendlier).
      proj = centred.dot(axis)
      axis *= -1.0 if proj.sum.negative?
      axis
    end

    # Bitmask over `rim_points2d` selecting points within +/- 120° of the
    # inferior pole, measuring theta from `si_dir2d` around `centroid`.
    def inferior_two_thirds_mask(rim_points2d, centroid:, si_dir2d:)
      n = rim_points2d.shape[0]
      mask = Numo::Bit.zeros(n)
      ix = si_dir2d[0]
      iy = si_dir2d[1]
      # Perpendicular axis (for atan2 in the SI-anchored frame).
      px = -iy
      py = ix
      n.times do |i|
        dx = rim_points2d[i, 0] - centroid[0]
        dy = rim_points2d[i, 1] - centroid[1]
        along = dx * ix + dy * iy   # +1 at inferior pole
        perp  = dx * px + dy * py
        theta = Math.atan2(perp, along)
        mask[i] = 1 if theta.abs <= PICO_HALF_ANGLE_RAD
      end
      mask
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

    # Walk the reference circle in the rim-plane 2D frame. For each angular
    # bin, cast a ray FROM the centre OUTWARD past the rim, and record the
    # maximum radius at which bone is still present.
    #
    # When `pico_sector` is nil, bone-loss is the full-circle area integral
    #   loss = (1/bins) * sum_b (1 - (max_bone_r_b / r)^2)
    # which is the polar-Jacobian form of (circle_area - bone_inside) / circle_area.
    #
    # When `pico_sector = { si_dir2d:, half_angle: }` is given, the integral
    # sums ONLY over bins inside +/- half_angle of the inferior pole (so the
    # natural pear-narrowing on the superior 1/3 doesn't get scored as loss),
    # and the denominator is the count of inferior-sector bins. This is the
    # clinical Pico method: bone loss is the missing area inside the
    # inferior-2/3 sector of the reference circle.
    def walk_circle_for_defect(circle:, projector:, affine:, scapula_mask:,
                               bins:, tolerance_mm:, pico_sector: nil)
      occupancy = Numo::Bit.zeros(bins)
      max_bone_r_arr = Array.new(bins, -1.0)
      affine_n = Affine.to_narray(affine)
      inv = Numo::Linalg.inv(affine_n)
      shape = scapula_mask.shape

      voxel_mm = Affine.voxel_spacing_mm(affine_n)
      step = [voxel_mm * 0.5, 0.25].max
      r_max = circle.radius + tolerance_mm
      n_samples = (r_max / step).ceil + 1

      normal = projector.normal
      n_through = 5
      through_spacing = voxel_mm

      # If we have a Pico sector, precompute which bins are "inferior".
      inferior_bin = nil
      if pico_sector
        inferior_bin = Numo::Bit.zeros(bins)
        ix = pico_sector[:si_dir2d][0]
        iy = pico_sector[:si_dir2d][1]
        half = pico_sector[:half_angle]
        bins.times do |b|
          theta_bin = (b.to_f / bins) * 2.0 * Math::PI
          # Unit vector at this bin angle (relative to circle centre, in (u,v)).
          bx = Math.cos(theta_bin)
          by = Math.sin(theta_bin)
          # Angle of this bin relative to the inferior pole.
          along = bx * ix + by * iy
          perp  = bx * (-iy) + by * ix
          phi = Math.atan2(perp, along)
          inferior_bin[b] = 1 if phi.abs <= half
        end
      end

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
        max_bone_r_arr[b] = max_bone_r
        occupancy[b] = 1 if max_bone_r >= circle.radius - tolerance_mm
      end

      # Area-based loss fraction.
      r = circle.radius
      area_total = Math::PI * (r ** 2)
      loss_sum = 0.0
      denom = 0
      max_bone_r_arr.each_with_index do |mbr, b|
        if pico_sector
          next unless inferior_bin[b] == 1
        end
        eff = (mbr.negative? ? 0.0 : [mbr, r].min)
        loss_sum += 1.0 - (eff / r) ** 2
        denom += 1
      end
      denom = bins if denom.zero?
      fraction = loss_sum / denom
      fraction = 0.0 if fraction.negative?
      fraction = 1.0 if fraction > 1.0
      # For the Pico method, area_lost is the lost fraction of the inferior-2/3
      # sector area. For the full-circle method it's the lost fraction of the
      # whole disc. Either way `bone_loss_fraction` is dimensionless and
      # `area_total` is the denominator's geometric area.
      total_for_area = pico_sector ? (Math::PI * r ** 2 * (denom.to_f / bins)) : area_total
      area_lost = fraction * total_for_area

      gap = DefectArc.largest_gap(occupancy, bins)
      angular_span = gap[:span] * (2.0 * Math::PI / bins)
      arc_length = angular_span * r

      start_angle = gap[:start] * (2.0 * Math::PI / bins)
      end_angle   = gap[:end_excl] * (2.0 * Math::PI / bins)
      DefectArc::Arc.new(start_angle, end_angle, arc_length, fraction,
                         area_lost, total_for_area, angular_span, occupancy)
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

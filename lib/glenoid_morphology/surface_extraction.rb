# frozen_string_literal: true

require "numo/narray"

module GlenoidMorphology
  # Extract the surface (boundary) voxels of a 3D boolean mask.
  #
  # A voxel is a surface voxel iff it is `true` AND at least one of its
  # six face-neighbours is `false` (or outside the array). We avoid building
  # a full erosion kernel by shifting the mask along each axis with `false`
  # padding and OR-ing the "neighbour is false" tests.
  module SurfaceExtraction
    module_function

    # @param mask [Numo::Bit] 3D boolean array (shape [D, H, W])
    # @return [Numo::Int64] array of shape [N, 3] of (i, j, k) voxel indices
    def surface_voxels(mask)
      raise ArgumentError, "mask must be 3-D" unless mask.ndim == 3

      bit = mask.cast_to(Numo::Bit)
      bool = bit.eq(1) # safe Numo::Bit

      boundary = Numo::Bit.zeros(*bit.shape)

      # For each of the six face-neighbours, OR in a "neighbour is background" map,
      # then AND with the mask itself.
      [
        shift(bool, axis: 0, by: 1),
        shift(bool, axis: 0, by: -1),
        shift(bool, axis: 1, by: 1),
        shift(bool, axis: 1, by: -1),
        shift(bool, axis: 2, by: 1),
        shift(bool, axis: 2, by: -1)
      ].each do |neighbour|
        boundary = boundary | (bool & neighbour.eq(0))
      end

      indices_of_true(boundary)
    end

    # Shift a Numo::Bit array along one axis, padding the vacated slab with 0.
    # by = +1 means "what was at index i is now at index i+1"; i.e. we look at
    # i-1 to ask "is my neighbour on the -axis side?". Conceptually we want
    # "the array of my +axis neighbours", which is shift(by: -1).
    def shift(arr, axis:, by:)
      out = Numo::Bit.zeros(*arr.shape)
      src_slicer = [true, true, true]
      dst_slicer = [true, true, true]
      n = arr.shape[axis]
      if by.positive?
        src_slicer[axis] = 0...(n - by)
        dst_slicer[axis] = by...n
      else
        b = -by
        src_slicer[axis] = b...n
        dst_slicer[axis] = 0...(n - b)
      end
      out[*dst_slicer] = arr[*src_slicer]
      out
    end

    # Return the [N, 3] integer index array of every true voxel in a Numo::Bit.
    def indices_of_true(bit)
      flat_idx = bit.where # Numo::Int32 of flat indices
      n = flat_idx.size
      return Numo::Int64.zeros(0, 3) if n.zero?

      _d, h, w = bit.shape
      hw = h * w
      i = flat_idx / hw
      r = flat_idx % hw
      j = r / w
      k = r % w
      out = Numo::Int64.zeros(n, 3)
      out[true, 0] = i
      out[true, 1] = j
      out[true, 2] = k
      out
    end

    # Given an [N, 3] voxel index array and a direction vector in voxel-space,
    # keep only the points whose outward-pointing position relative to the
    # centroid has a positive dot product with `direction`. Used to retain
    # the glenoid-facing surface (the side facing the humerus).
    #
    # @param voxels  [Numo::Int64] [N, 3]
    # @param direction [Array<Float>] length 3 in voxel space
    # @return [Numo::Int64] filtered subset
    def filter_facing(voxels, direction:)
      n = voxels.shape[0]
      return voxels if n.zero?

      v = Numo::DFloat.cast(voxels)
      centroid = v.mean(axis: 0)
      rel = v - centroid.reshape(1, 3)
      d = Numo::DFloat.cast(direction)
      d /= Math.sqrt((d ** 2).sum)
      dots = rel.dot(d)
      keep = dots > 0
      idx = keep.where
      voxels[idx, true]
    end

    # Largest 6-connected component of a boolean mask, returned as a Numo::Bit
    # of the same shape with only the component voxels set. Implemented as an
    # iterative stack-based flood fill in pure Ruby — slow-ish (a few seconds
    # for ~500k voxels) but dependency-free.
    #
    # Real segmentation masks routinely have spurious blobs (an acromion tip,
    # a stray rib voxel cluster, etc.) that hijack the PCA plane fit. Keeping
    # only the largest connected component drops them before any geometry.
    def largest_component(mask)
      bit = mask.cast_to(Numo::Bit)
      shape = bit.shape
      raise ArgumentError, "mask must be 3-D" unless shape.length == 3

      flat_true = bit.where.to_a
      return Numo::Bit.zeros(*shape) if flat_true.empty?

      in_set = {}
      flat_true.each { |fi| in_set[fi] = true }
      visited = {}
      hw = shape[1] * shape[2]
      w  = shape[2]

      best_size = 0
      best_set  = nil
      flat_true.each do |start|
        next if visited[start]

        stack = [start]
        visited[start] = true
        component = []
        while (cur = stack.pop)
          component << cur
          i = cur / hw
          r = cur % hw
          j = r / w
          k = r % w
          [[i + 1, j, k], [i - 1, j, k],
           [i, j + 1, k], [i, j - 1, k],
           [i, j, k + 1], [i, j, k - 1]].each do |ni, nj, nk|
            next if ni < 0 || nj < 0 || nk < 0
            next if ni >= shape[0] || nj >= shape[1] || nk >= shape[2]

            fi = ni * hw + nj * w + nk
            next if visited[fi]
            next unless in_set[fi]

            visited[fi] = true
            stack << fi
          end
        end
        if component.size > best_size
          best_size = component.size
          best_set  = component
        end
      end

      out = Numo::Bit.zeros(*shape)
      best_set.each do |fi|
        i = fi / hw
        rem = fi % hw
        j = rem / w
        k = rem % w
        out[i, j, k] = 1
      end
      out
    end

    # Keep only voxels (in [N, 3] index space) whose physical distance to
    # `anchor_mm` is at most `min_distance_mm + window_mm`, where
    # `min_distance_mm` is the smallest distance from any voxel in the set to
    # the anchor. Used to crop a whole-scapula surface to just the glenoid
    # surface region by anchoring on the humerus centroid.
    #
    # @param voxels [Numo::Int64] [N, 3] voxel indices
    # @param voxels_mm [Numo::DFloat] [N, 3] same voxels in physical mm
    # @param anchor_mm [Numo::DFloat] [3] anchor point in physical mm
    # @param window_mm [Float] keep voxels with d <= d_min + window_mm
    # @return [Numo::Int64] filtered subset of `voxels`
    def filter_near_anchor(voxels, voxels_mm:, anchor_mm:, window_mm:)
      n = voxels.shape[0]
      return voxels if n.zero?

      ax, ay, az = anchor_mm[0], anchor_mm[1], anchor_mm[2]
      dx = voxels_mm[true, 0] - ax
      dy = voxels_mm[true, 1] - ay
      dz = voxels_mm[true, 2] - az
      dist = Numo::NMath.sqrt(dx * dx + dy * dy + dz * dz)
      d_min = dist.min
      keep = dist <= (d_min + window_mm)
      idx = keep.where
      voxels[idx, true]
    end
  end
end

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
  end
end

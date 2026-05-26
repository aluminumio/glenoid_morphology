# frozen_string_literal: true

require "numo/narray"

module GlenoidMorphology
  # Voxel-index <-> physical-mm transforms backed by a 4x4 NIfTI-style affine.
  #
  # The affine maps voxel index [i, j, k, 1] to physical [x, y, z, 1] in mm.
  # We accept either a Numo::DFloat [4,4] or a plain Ruby [[..]] array.
  module Affine
    module_function

    def to_narray(affine)
      return affine if affine.is_a?(Numo::DFloat) && affine.shape == [4, 4]

      Numo::DFloat.cast(affine).tap do |a|
        raise ArgumentError, "affine must be 4x4" unless a.shape == [4, 4]
      end
    end

    # Apply affine to an [N, 3] Numo array of voxel indices, returning [N, 3] mm.
    def voxels_to_mm(voxels, affine)
      affine = to_narray(affine)
      v = Numo::DFloat.cast(voxels)
      n = v.shape[0]
      hom = Numo::DFloat.ones(n, 4)
      hom[true, 0..2] = v
      out = hom.dot(affine.transpose)
      out[true, 0..2]
    end

    # Approximate isotropic voxel spacing (mean of column norms of the rotation part).
    def voxel_spacing_mm(affine)
      a = to_narray(affine)
      cols = [0, 1, 2].map { |i| Math.sqrt((a[0..2, i] ** 2).sum) }
      cols.sum / 3.0
    end

    # Identity affine with given isotropic spacing (handy for tests).
    def identity(spacing_mm: 1.0)
      a = Numo::DFloat.eye(4)
      a[0, 0] = spacing_mm
      a[1, 1] = spacing_mm
      a[2, 2] = spacing_mm
      a
    end
  end
end

# frozen_string_literal: true

require_relative "glenoid_morphology/version"
require_relative "glenoid_morphology/measure"

# Geometric glenoid bone-loss measurement from a 3D scapula label volume.
#
# The pipeline (see README and doc/method.md) implements the en-face /
# best-fit-circle method described in arxiv 2511.14083:
#
#   1. Extract the lateral, glenoid-facing surface of the scapula mask.
#   2. Fit an en-face plane to those surface voxels via PCA.
#   3. Project the surface onto the plane → 2D rim candidates.
#   4. RANSAC + least-squares-fit a circle to the rim (the intact glenoid).
#   5. Identify the largest empty arc → the bone defect.
#   6. Compute bone-loss percent from the circular segment area.
#
# All math runs on Numo::NArray + Numo::Linalg; no neural networks here.
module GlenoidMorphology
  class Error < StandardError; end

  # @param scapula_mask [Numo::Bit] 3D boolean mask (the scapula from `shoulder_segmenter`)
  # @param affine       [Array<Array<Float>>, Numo::DFloat] 4x4 voxel → mm transform
  # @param humerus_mask [Numo::Bit, nil] optional, used to orient the glenoid-facing side
  # @param side         [:left, :right, nil] anatomical side, inferred if omitted
  # @param opts         [Hash] tuning knobs (see Measure::DEFAULT_OPTS)
  # @return [Measurement]
  def self.measure(scapula_mask:, affine:, **opts)
    Measure.call(scapula_mask: scapula_mask, affine: affine, **opts)
  end
end

# glenoid_morphology

Geometric **glenoid bone-loss measurement** from shoulder CT segmentations, in
pure Ruby. Takes the per-voxel scapula label produced by
[`shoulder_segmenter`](https://github.com/aluminumio/shoulder_segmenter) and
returns:

- the en-face plane of the glenoid,
- the best-fit circle approximating the intact glenoid rim,
- the defect arc (start angle, end angle, arc length), and
- bone-loss expressed as a percent of the original rim area.

This is a Ruby port of the geometry pipeline described in
[arxiv 2511.14083](https://arxiv.org/abs/2511.14083) (en-face plane / best-fit
circle / chord-segment area). It is **pure math** — RANSAC, PCA via
`numo-linalg.eigh`, Kasa-method circle fitting, mask ray-casting — no neural
networks.

## Public API

```ruby
require "glenoid_morphology"

result = GlenoidMorphology.measure(
  scapula_mask: numo_bit_3d,           # Numo::Bit shape [D, H, W]
  humerus_mask: numo_bit_3d_optional,  # optional, used to orient the glenoid-facing side AND crop the surface to the glenoid region
  affine:       voxel_to_mm_4x4,       # 4x4 voxel→mm transform (from the NIfTI header)
  side:         :right,                # :left or :right; inferred if omitted
  glenoid_window_mm: 30.0,             # crop radius around humerus when humerus_mask given (default 30 mm)
  largest_component: true              # drop disconnected blobs in the scapula mask before extracting surface (default true)
)

result.bone_loss_percent   # => e.g. 18.4
result.fitted_circle       # => { center: [x,y,z], radius_mm: 23.8, normal: [...] }
result.en_face_plane       # => { a:, b:, c:, d: }  (ax + by + cz = d in mm)
result.defect_arc          # => { start_angle:, end_angle:, arc_length_mm: }
result.confidence          # => 0.0..1.0 (low if fit residual is large)
```

See `doc/method.md` for the full algorithm.

## Installation

```ruby
# Gemfile
gem "glenoid_morphology"
```

### `numo-linalg` on macOS arm64

`numo-linalg` needs a LAPACK backend. On Apple silicon (M1/M2/M3/M4) the
clean install path is:

```sh
brew install openblas
bundle config build.numo-linalg --with-openblas-dir=/opt/homebrew/opt/openblas
bundle install
```

Without that flag the gem will attempt to find Accelerate/MKL and may fail
to link. This recipe was verified on macOS 14 / arm64 against
`openblas 0.3.33` and `numo-linalg 0.1.7`.

## Algorithm

Implements the geometry pipeline of [arxiv:2511.14083][paper] without any
neural-network involvement:

1. **Isolate the glenoid-facing surface** of the scapula mask via
   morphological-gradient surface extraction + directional filter.
2. **Fit the en-face plane** by PCA — the eigenvector of the surface-point
   covariance matrix with the smallest eigenvalue is the plane normal.
3. **Project the surface** onto the plane → 2D rim candidates.
4. **Extract the rim** as the 2D convex hull of the projection.
5. **RANSAC + LSQ best-fit circle** through the rim points
   (Andrew's chain hull + Kasa algebraic LSQ with iterative reweighting).
6. **Walk the fitted circle** and check, at each angle, whether bone
   exists near the rim radius. Contiguous "no-bone" runs are the defect arc.
7. **Bone loss** = circular segment area / disc area
   = `(theta - sin(theta)) / (2 * pi)`.

[paper]: https://arxiv.org/abs/2511.14083

## Where it sits in the imaging stack

```
       DICOM series
            ↓  (dicom_seg-ruby)
       NIfTI volume
            ↓  (shoulder_segmenter)
   per-voxel bone labels
            ↓  (glenoid_morphology — THIS GEM)
   bone-loss percent + en-face plane + best-fit circle
```

This gem only depends on `numo-narray` and `numo-linalg`. It deliberately
does **not** depend on torch-rb, nifti-ruby, or Rails — it's a self-contained
geometry library, runnable with a synthetic mask and no external assets.

## Synthetic specs

The Tier-1 unit specs cover every math primitive (plane fit, circle fit,
defect arc, surface extraction, projector, affine). The Tier-2 end-to-end
specs build synthetic voxelised glenoid masks with known chord defects and
verify the algorithm recovers the bone-loss percent within a few percent
of the geometric expected value. See `spec/support/synthetic.rb`.

At 1 mm voxel spacing, the practical detection floor for a ~24 mm radius
synthetic glenoid is ~4 percentage-points of bone loss. Real CT data with
sub-millimetre spacing will be sharper.

## License

MIT.

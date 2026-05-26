# Method

This document describes the geometric pipeline implemented in
`GlenoidMorphology.measure`. It is a Ruby port of the en-face / best-fit
circle method described in [arxiv 2511.14083][paper] for measuring anterior
glenoid bone loss from a 3D scapula segmentation.

[paper]: https://arxiv.org/abs/2511.14083

## Inputs

- **Scapula mask**: a 3D boolean volume (`Numo::Bit` shape `[D, H, W]`)
  identifying scapula voxels. Typically this is one channel of the
  `shoulder_segmenter` output.
- **Affine** (`4x4` `Numo::DFloat`): the voxel-index → physical-mm
  transform from the underlying NIfTI volume.
- **Humerus mask** (optional, recommended): used to disambiguate which
  face of the scapula is the glenoid (the lateral, articular surface
  facing the humerus). When provided, the glenoid-facing direction is
  `centroid(humerus) - centroid(scapula)` in voxel space.

## Pipeline

### 1. Surface extraction

A voxel is a surface voxel iff it is `true` in the mask AND has at least
one face-neighbour that is `false` (or outside the array). Equivalent to
`mask & !erode(mask)` for a 6-connected erosion kernel, but computed by
shifting the mask along each axis and OR-ing the boundary masks.

### 2. Glenoid-facing filter

The output of step 1 is the *whole* scapula surface; we only want the
lateral, articular face. We keep only points whose displacement from the
scapula-surface centroid has a positive dot product with the
glenoid-facing direction.

### 3. En-face plane fit (PCA)

Compute the `3x3` covariance matrix of the (centred) surface points and
take `Numo::Linalg.eigh`. The eigenvector with the smallest eigenvalue is
the plane normal; the plane passes through the centroid.

Returned as `(a, b, c, d)` with `(a, b, c)` a unit vector and `a*x + b*y +
c*z = d` in mm.

### 4. 2D projection

Pick an orthonormal in-plane basis `(u, v)` deterministically:

- `u` = unit projection of the world x-axis onto the plane (falling back to
  the world y-axis if x is nearly parallel to the normal),
- `v` = `normal × u`.

Project each 3D point onto `(u, v)` to obtain its 2D rim-plane coordinates.

### 5. Rim extraction (convex hull)

The articular surface of an intact OR chord-defect glenoid is convex in
the en-face plane, so the rim is exactly the **convex hull** of the
projected point cloud. We use Andrew's monotone-chain in
`O(n log n)`. Per-angle max-radius binning was tried first but produces
bimodal radius distributions for chord-cut discs that bias the circle fit;
the hull is cleaner.

### 6. Best-fit circle (RANSAC + Kasa LSQ + IRLS)

We do three passes:

1. **RANSAC** on the top-radius quantile of the hull, seeded with three
   exact-circle samples and `Kasa` LSQ refinement on the inlier set.
2. **RANSAC** on all hull points, with a wider iteration budget.
3. Keep whichever of (1, 2) has the larger radius — the chord side biases
   the radius downward, so the larger-radius fit is the one less corrupted
   by the defect.
4. **Iteratively reweighted Kasa** — re-fit using only points whose signed
   radial residual is `>= -tolerance`, repeat until centre+radius converge.
   This pulls the fit fully onto the intact arc.

Kasa's algebraic LSQ solves `min sum_i (x_i^2 + y_i^2 - 2*a*x_i - 2*b*y_i -
c)^2` where the circle is `(x-a)^2 + (y-b)^2 = r^2` with
`r^2 = a^2 + b^2 + c`. Closed-form, fast, well-conditioned for points that
span more than ~90° of arc.

### 7. Defect-arc detection by mask ray-casting

For each of `defect_bins` angular bins around the fitted circle:

- Cast a ray from the fitted centre outward in the rim plane.
- At each radial sample, lift back into 3D mm, transform to voxel index
  via the inverse affine, and check the mask.
- Sample a few offsets along the plane normal to handle small drift between
  the fitted plane and the true disc plane.
- Track the maximum radius at which bone is still present.
- The bin is intact iff `max_bone_radius >= fitted_radius - tolerance`.

The largest contiguous run of "no-bone" bins is the principal defect arc.

### 8. Bone-loss percent

The lost area is a polar-coordinate area integral of "missing bone inside
the fitted circle":

```
bone_loss_fraction = (1 / N_bins) * sum_b (1 - (max_bone_r_b / r)^2)
```

where `max_bone_r_b` is the largest radius at which bone is present along
the ray at angular bin `b`, and `r` is the fitted circle radius. This is
the polar-Jacobian form of `(circle_area - bone_area_inside_circle) / circle_area`.

For a clean chord defect (synthetic case) this reduces to the standard
circular-segment formula `(theta - sin(theta)) / (2*pi)` where `theta` is
the angular span of the chord. For a real glenoid the integral is also
robust to small mismatches between the fitted circle and the (mildly
elliptical) glenoid outline.

## Tunable knobs (`GlenoidMorphology::Measure::DEFAULT_OPTS`)

| key                          | default | meaning                                                                                |
|------------------------------|---------|----------------------------------------------------------------------------------------|
| `method`                     | `:pico` | `:pico` (inferior-2/3 reference circle, Baudi 2005) or `:full_circle` (legacy)         |
| `humerus_mask`               | nil     | optional Numo::Bit; enables glenoid-facing direction inference + proximity surface crop |
| `glenoid_window_mm`          | 30.0    | when `humerus_mask` is given, keep only surface voxels within `d_min + window_mm` of the humerus centroid |
| `largest_component`          | true    | strip disconnected blobs in the scapula mask before extracting surface                |
| `ransac_iterations`          | 200     | RANSAC iterations for circle fit                                                       |
| `ransac_threshold_mm`        | 1.5     | mm tolerance for RANSAC inliers                                                        |
| `rim_radial_tolerance_mm`    | 1.0     | mm tolerance for "on-rim" occupancy bit in step 7                                      |
| `defect_bins`                | 360     | angular bin count (1° per bin)                                                         |
| `seed`                       | nil     | RNG seed for reproducible RANSAC                                                       |

## Confidence

A heuristic blend of (a) RANSAC inlier ratio, (b) Kasa residual relative
to fitted radius, (c) number of rim points, and (d) plane fit residual.
Returned in `[0, 1]`; values below ~0.5 indicate a poor fit or insufficient
data and the bone-loss percent should not be trusted.

## Numerical notes

- All math is `Numo::DFloat` (double precision).
- `atan2` is not vectorised in Numo, so angular calculations are done in a
  Ruby loop. For typical scapula sizes (~few thousand rim voxels) this is
  fine; the heavy lifting (PCA, RANSAC, mask sampling) is the bottleneck.
- The convex hull is implemented on plain Ruby tuples (small N).

## Pico inferior-2/3 method (v0.3.0 default)

The healthy glenoid is pear-shaped (~30 mm superior-inferior × ~18 mm
anterior-posterior). The lower 2/3 of the rim is circular; the superior
1/3 narrows. Fitting a single circle to the WHOLE rim averages the SI
and AP dimensions and then the area-integral bone-loss metric flags the
ellipse-vs-circle shape difference as bone loss (~55% on a healthy
glenoid). The Pico method fits the reference circle to only the
inferior-2/3 of the rim — which is genuinely circular — then measures
bone loss only over that same 240° sector.

Citation: Baudi P, Righi P, Bolognesi D, Campochiaro G, Rebuzzi M,
Matino G, Catani F. *How to identify and calculate glenoid bone deficit*.
Chir Organi Mov. 2005; 90(2):145-152. (Variants are also published as
"Pico" / "best-fit inferior circle" in Sugaya 2003 and subsequent
literature.)

Pipeline additions on top of steps 1-6 above:

1. Compute the inferior direction in the en-face plane's 2D basis. The
   world-space inferior unit vector is `-superior`, derived from the
   NIfTI affine: pick the affine column most aligned with world +Z (the
   image axis closest to anatomical superior), negate it. For a standard
   RAS+ affine this is `[0, 0, -1]`.
2. Project the 3D inferior unit onto the en-face plane (subtract the
   normal component) and re-normalise; express in `(u, v)` basis. That
   gives a 2D unit vector pointing inferior in projected coordinates.
3. Filter rim points to those within `±120°` of the inferior pole
   (measuring from the projected inferior direction around the rim
   centroid). Sanity check: if fewer than 20% of rim points survive,
   raise; the SI direction was probably wrong.
4. Fit the **Pico reference circle** to just those inferior-2/3 points
   using the same RANSAC + Kasa + IRLS pass as the full-rim circle.
5. Walk the full 360° of the Pico circle as before, but in the
   area-integral sum at the end of step 7 above, sum and divide only
   over the inferior-2/3 bins. The natural pear-narrowing on the
   superior third is then NOT scored as bone loss.

### Edge case: SI direction is nearly parallel to the plane normal

For an unusual scanner orientation (glenoid facing straight up/down),
the projection of world-inferior onto the en-face plane is degenerate.
In that case we fall back to the **longest principal axis of the 2D rim
point cloud** as the SI proxy, sign-disambiguated to point toward the
majority of the rim points (i.e. away from the narrowed superior tip).
A warning is emitted.

### `method: :full_circle` (legacy)

The pre-Pico algorithm is still available via the `method: :full_circle`
kwarg. The struct field `legacy_bone_loss_percent` always carries the
full-circle bone-loss for back-comparison, regardless of which method was
selected as the primary.

## Known issue: real-mask bone-loss is sensitive to seg quality (v0.3.0)

On the `shoulder_segmenter` `scapula_left` mask we instrumented with
`script/debug_real_mask.rb` (`tmp/segmentation_run_appendicular/labels.nii.gz`),
the Pico method drops bone-loss from the full-circle ~60% but settles
around 40-50% rather than the expected ~0-10% for a healthy shoulder.

Drilling in: the SI direction projection from the affine is correct
(`world_inferior ≈ [0.1, 0.13, -0.99]`, projection onto the en-face
plane is non-degenerate, no fallback fires), and the rim point cloud has
the expected ~30mm × 18mm pear-shaped extent in projected coordinates.
The problem is upstream: the convex-hull rim has only ~25 points and
substantial gaps on one side of the inferior arc. With a small subset of
points, the inferior-2/3 fit collapses onto a tight curve and reports
`r ≈ 6-8 mm` instead of the clinical ~12-13 mm. Sweeping the proximity
crop window pushes the radius back to ~12 mm but at the cost of pulling
non-glenoid scapular surface into the rim.

Likely follow-ups:
- Radial smoothing or morphological closing on the mask before surface
  extraction.
- Denser rim sampling (per-angle max-radius after the convex hull, with
  outlier rejection on the chord side).
- Cross-validate the segmentation against a CT viewer; the mask may
  genuinely be truncated.

This is upstream of the Pico geometry; the Pico method works correctly
on the synthetic pear glenoid (`spec/pico_method_spec.rb`).

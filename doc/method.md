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

## Limitations vs. the paper

- We use Kasa for circle fitting (the paper notes Pratt as an alternative
  for sub-degree-arc inputs). Kasa is sufficient for the wide arcs we see
  in real glenoids (rim covers >270°).
- We assume the defect is a single contiguous chord-segment removal. The
  algorithm reports only the LARGEST defect arc.
- Real glenoid measurements may also use the inferior-circle method (paired
  inferior-quadrant circle vs. full-circle). That isn't implemented yet —
  if needed it's a small addition built on the same circle-fit primitives.

## Known issue: bone-loss over-estimation on real masks (v0.2.0)

On a real `scapula_left` mask from `shoulder_segmenter`, the radius now
fits the glenoid (~13 mm) but `bone_loss_percent` still reports ~55% for
an apparently healthy glenoid. Root cause: the real glenoid is pear/
elliptical-shaped (~30 mm SI × ~18 mm AP), and the best-fit circle picks
a radius between those dimensions, so the area-integral metric flags
shape-vs-circle mismatch as bone-loss. The fix here is the Pico
inferior-2/3 circle method (fit only the inferior arc) — deferred to a
follow-up.

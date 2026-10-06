# camcal236 — planar-target camera calibration SDK

C++20 camera calibration library (namespace `camcal236`) built from scratch on
Eigen 3.4.0 (headers in `third_party/eigen3`). No frontend, no HTTP, no
corner detection, no distortion estimation: only the ideal pinhole model (or
pre-undistorted observations) with zero skew.

## Build

```sh
make            # builds libcamcal236.a and the demo binary
make test       # builds and runs the known-parameter multi-view example
```

Requires g++ 11.4 (`-std=c++20`) and `-Ithird_party/eigen3`.

## API (include/camcal236/camcal.hpp)

- `calibrate(target_points, views, options)` — full calibration.
- `project(fx, fy, cx, cy, R, t, point_mm)` — projection for the same model.

Input: up to 48 planar target points (millimetres, z = 0 plane, x right /
y down, unique finite IDs) and 3–6 views. Observations are associated to
target points by ID; missing points are allowed, but every view needs at
least 8 observations. Image width/height must be positive integers and every
observation must lie inside the image. Unknown IDs, duplicate observations,
non-finite values and (near-)collinear geometry are rejected; inputs are
never modified.

Output (`CalibrationResult`): shared intrinsics `fx, fy, cx, cy` (skew = 0),
per-view target-to-camera rotation (projected to SO(3)) and translation,
per-observation prediction and residual, per-view and overall RMS
(root-mean-square of per-point squared 2D distances, pixels). All outputs are
recomputed from the final parameters.

## Method

1. **Validation** (`src/validation.cpp`) — all input checks above.
2. **Homography** (`src/homography.cpp`) — per-view normalized (Hartley) DLT.
3. **Geometric initialization** (`src/geometry_init.cpp`) — Zhang-style
   closed-form intrinsics from the orthogonality constraints of all views,
   then per-view extrinsics. Rotations are projected to SO(3) via SVD and the
   sign is chosen so observed target points have positive depth. If the
   intrinsics are not identifiable (e.g. near-parallel target planes, no
   positive focal solution), calibration fails explicitly and delivers no
   parameters — no fixed-focal fallback is used.
4. **Joint optimization** (`src/optimizer.cpp`) — Levenberg–Marquardt over
   the shared intrinsics and all poses jointly, minimizing the sum of squared
   2D pixel residuals over every observation (no per-view fit-and-average).
   `fx, fy` are parameterized in log-space to stay positive; rotations are
   updated on SO(3) via the exponential map. Only objective-non-increasing
   updates are accepted. Iteration cap and convergence tolerance are
   configurable via `SolverOptions`. If the cap is reached, the last valid
   estimate is returned with `converged = false`; if initialization fails,
   no parameters are delivered (`success = false`).

## Degeneracy tolerances

- `collinearity_tolerance` (default `1e-6`): a point set is rejected as
  (near-)collinear when the ratio of the smaller to the larger singular value
  of its centered 2D coordinates falls below this (scale-invariant). Applied
  to the full target and to the subset observed in each view.
- Homography degeneracy: the 8th singular value of the DLT system must exceed
  `1e-12`.
- Intrinsic identifiability: the null space of the Zhang constraint matrix
  must be 1-dimensional — the second-smallest singular value must exceed
  `1e-10` times the largest, otherwise the views are too similar
  (near-parallel planes) and calibration fails.
- `depth_epsilon` (default `1e-9` mm): minimum mean target depth accepted
  when disambiguating the extrinsic sign.

## Scope and limitations

- Ideal pinhole only; observations must be distortion-free or pre-undistorted.
- No corner detection, no distortion estimation, no image I/O.
- Planar single-target calibration with shared intrinsics and zero skew.
- 3–6 views, up to 48 target points, at least 8 observations per view.

## Example

`examples/demo.cpp` synthesizes a 8x6 grid (25 mm spacing, 48 points),
renders 4 noisy views (0.15 px noise, ~10% of points randomly missing per
view) from known intrinsics/poses, calibrates, and prints the comparison.
Typical result: intrinsics within ~0.2% of ground truth, overall RMS at the
noise level.

## License

Eigen is licensed under MPL2 (see `third_party/eigen3/COPYRIGHT`).


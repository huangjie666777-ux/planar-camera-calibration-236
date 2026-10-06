# camcal236 — planar-target camera calibration SDK

C++20 / Eigen 3.4.0, no frontend, no HTTP. Namespace camcal236.

## Build and test

    make            # builds build/libcamcal236.a, demo and selftest
    make test       # runs tests/selftest.cpp and examples/multiview_demo.cpp

Compile flags: -std=c++20 -Ithird_party/eigen3 -Iinclude, link build/libcamcal236.a.

## API

See include/camcal236/types.hpp and include/camcal236/calibration.hpp.

- Calibrate(target, views, image, options) — full pipeline: validation,
  geometric initialisation, joint refinement. Never modifies its inputs.
- Project(intrinsics, pose, target_xy_mm) — pinhole projection of a
  target-plane point (mm, z = 0) into pixels.

Input contract:

- Up to 48 target points (mm, x right / y down on the plane), unique ids.
- 3-6 views, at least 8 observations per view, associated by target id;
  missing points per view are allowed.
- Image width/height are positive integers; observations must lie inside
  0 <= u < width, 0 <= v < height and be finite.
- Rejected: unknown ids, duplicate target ids, duplicate observations per
  view, non-finite values, collinear geometry.
- Ideal pinhole only (zero skew, no distortion): observations must be
  ideal or pre-undistorted. No corner detection, no distortion estimation.

Output (Output): on validation/initialisation failure ok == false and no
parameters are delivered. Otherwise result carries shared intrinsics
(fx, fy, cx, cy), per-view target-to-camera poses (R on SO(3), t), per-point
predictions and residuals, per-view RMS and total RMS — all recomputed from
the final parameters, RMS = sqrt(mean of squared 2D distances). If the
optimiser hits the iteration cap, the last valid estimate is kept and
converged == false.

## Method

1. Validation (src/validation.cpp).
2. Homographies (src/homography.cpp): Hartley-normalised DLT per view.
3. Initialisation (src/initialization.cpp): Zhang-style orthogonality
   constraints on B = K^-T K^-1 stacked over all views; closed-form
   intrinsics (skew dropped, model fixes it to zero); per-view poses from
   K^-1 H, projected onto SO(3) via SVD, with the sign ambiguity resolved by
   requiring positive target depth. If the constraint system is rank-deficient
   or the focal lengths come out non-positive, initialisation fails explicitly
   ("intrinsics not identifiable") — no fixed-focal-length fallback.
4. Joint refinement (src/optimizer.cpp): Levenberg-Marquardt over shared
   intrinsics and all view poses (Rodrigues + translation), minimising the sum
   of squared 2D pixel residuals of every observation of every view (no
   per-view fitting and averaging). Only cost-non-increasing updates are
   accepted; fx, fy are kept positive by step rejection; rotations stay on
   SO(3) by construction. Iteration cap and convergence tolerance are
   configurable via Options.

## Degeneracy tolerances

Defined in src/internal.hpp:

- kCollinearityTol = 1e-6: a point set (target points, or per-view observed
  pixels / target points) is rejected as collinear when the smallest singular
  value of its centred coordinate scatter is below tol^2 times the largest
  (per point).
- kSingularTol = 1e-12: the homography DLT system must have rank 8 and the
  intrinsic constraint system rank 5 (second-smallest singular value above
  this bound), otherwise the geometry is reported degenerate. The image of the
  absolute conic must be positive definite (sign test; its entries scale with
  1/f^2, so no absolute magnitude test is meaningful).

## Scope

Planar targets, 3-6 views, shared intrinsics, ideal pinhole. Not covered:
distortion estimation, corner detection, non-planar rigs, per-view intrinsics.

## Environment

C++20 is supported by GCC11.4.0. GNU Make4.3 is available. Eigen3.4.0 headers are in third_party/eigen3; compile with -std=c++20 -Ithird_party/eigen3. License is in third_party/eigen3/COPYRIGHT. Exact package provenance is in dependencies.lock.json.

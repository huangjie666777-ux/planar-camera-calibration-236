#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace camcal236 {

// A planar-target control point. Coordinates are millimetres on the target
// plane (x right, y down, z = 0). Ids must be unique within a target.
struct TargetPoint {
  int id = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;
};

// One pixel observation of a target point, associated by id.
struct Observation {
  int id = 0;
  double u = 0.0;  // pixel x (right)
  double v = 0.0;  // pixel y (down)
};

// Observations of one view. At least 8 observations per view are required.
struct View {
  std::vector<Observation> observations;
};

struct ImageSize {
  int width = 0;   // pixels, must be > 0
  int height = 0;  // pixels, must be > 0
};

// Shared pinhole intrinsics, skew fixed to zero. Only ideal pinhole or
// pre-undistorted observations are supported.
struct Intrinsics {
  double fx = 0.0;
  double fy = 0.0;
  double cx = 0.0;
  double cy = 0.0;
};

// Rigid transform mapping target-plane coordinates (z = 0) into the camera
// frame of one view: p_cam = R * [x, y, 0]^T + t.
struct Pose {
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  Eigen::Vector3d t = Eigen::Vector3d::Zero();
};

struct Options {
  int max_iterations = 200;      // LM iteration cap
  double convergence_tol = 1e-12;  // relative cost / step / gradient tolerance
};

struct PointResult {
  int id = 0;
  Eigen::Vector2d predicted = Eigen::Vector2d::Zero();  // pixels
  Eigen::Vector2d residual = Eigen::Vector2d::Zero();   // observed - predicted
};

struct ViewResult {
  Pose pose;
  double rms = 0.0;  // sqrt(mean squared 2D distance) over this view's points
  std::vector<PointResult> points;
};

struct CalibrationResult {
  bool converged = false;
  int iterations = 0;
  Intrinsics intrinsics;
  std::vector<ViewResult> views;
  double total_rms = 0.0;  // sqrt(mean squared 2D distance) over all points
};

struct Output {
  bool ok = false;
  std::string error;  // human-readable failure reason when ok == false
  CalibrationResult result;
};

}  // namespace camcal236

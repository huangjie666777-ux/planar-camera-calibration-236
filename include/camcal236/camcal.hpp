#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace camcal236 {

// A planar target point in millimetres, on the z = 0 target plane.
// x to the right, y down. IDs must be unique and finite.
struct TargetPoint {
    int id = 0;
    double x_mm = 0.0;
    double y_mm = 0.0;
    bool operator==(const TargetPoint&) const = default;
};

// One pixel observation of a target point, x right / y down, inside the image.
struct Observation {
    int id = 0;
    double u = 0.0;
    double v = 0.0;
    bool operator==(const Observation&) const = default;
};

// One view: positive integer image size plus its observations.
struct ViewInput {
    int width = 0;
    int height = 0;
    std::vector<Observation> observations;
};

struct SolverOptions {
    int max_iterations = 200;             // LM iteration cap
    double convergence_tolerance = 1e-10; // relative cost-decrease tolerance
    // Degeneracy tolerances (documented in README):
    double collinearity_tolerance = 1e-6; // spread of used target points vs scale
    double depth_epsilon = 1e-9;          // minimum accepted depth margin
};

struct ViewResult {
    Eigen::Matrix3d rotation;     // target -> camera
    Eigen::Vector3d translation;  // target -> camera, millimetres
    double rms = 0.0;             // per-view RMS of 2D residuals (pixels)
};

struct PointResult {
    int view_index = 0;
    int id = 0;
    Eigen::Vector2d predicted{0.0, 0.0};
    Eigen::Vector2d residual{0.0, 0.0};  // observed - predicted
};

struct CalibrationResult {
    bool success = false;     // false: no parameters delivered
    bool converged = false;   // false: iteration cap hit; last valid estimate kept
    std::string message;

    double fx = 0.0, fy = 0.0, cx = 0.0, cy = 0.0;  // skew fixed to zero
    std::vector<ViewResult> views;
    std::vector<PointResult> points;  // per observation, recomputed from final params
    double rms = 0.0;                 // overall RMS of 2D residuals (pixels)
};

// Full planar-target calibration. Does not modify its inputs.
CalibrationResult calibrate(const std::vector<TargetPoint>& target_points,
                            const std::vector<ViewInput>& views,
                            const SolverOptions& options = SolverOptions{});

// Projection interface for the same model (ideal pinhole, zero skew).
// point_mm lies on the target plane (z = 0 implied).
Eigen::Vector2d project(double fx, double fy, double cx, double cy,
                        const Eigen::Matrix3d& rotation,
                        const Eigen::Vector3d& translation,
                        const Eigen::Vector2d& point_mm);

}  // namespace camcal236

#pragma once

#include "camcal236/camcal.hpp"

#include <unordered_map>
#include <vector>

namespace camcal236::detail {

// Validated, ID-resolved working data. Built without modifying caller input.
struct Problem {
    std::vector<TargetPoint> points;                 // copied target points
    std::unordered_map<int, std::size_t> index_of;   // id -> index into points
    std::vector<std::vector<std::pair<std::size_t, Eigen::Vector2d>>> views;  // (point index, pixel)
    std::vector<ViewInput> raw_views;                // copied view metadata
};

// Returns empty string on success, else a human-readable rejection reason.
std::string validate(const std::vector<TargetPoint>& target_points,
                     const std::vector<ViewInput>& views,
                     const SolverOptions& options,
                     Problem& problem);

// Normalized-DLT homography from target plane (z=0) to pixels for one view.
// Returns false if the geometry is degenerate.
bool estimate_homography(const Problem& problem, std::size_t view,
                         Eigen::Matrix3d& H);

// Zhang-style multi-view intrinsics + per-view extrinsics initialization.
// Returns false with a message when intrinsics are not identifiable.
bool initialize(const Problem& problem, const SolverOptions& options,
                double& fx, double& fy, double& cx, double& cy,
                std::vector<Eigen::Matrix3d>& rotations,
                std::vector<Eigen::Vector3d>& translations,
                std::string& error);

// Joint Levenberg-Marquardt refinement of shared intrinsics and all poses.
// Only accepts non-increasing updates. Sets converged=false if the cap is hit.
void refine(const Problem& problem, const SolverOptions& options,
            double& fx, double& fy, double& cx, double& cy,
            std::vector<Eigen::Matrix3d>& rotations,
            std::vector<Eigen::Vector3d>& translations,
            bool& converged, int& iterations);

}  // namespace camcal236::detail


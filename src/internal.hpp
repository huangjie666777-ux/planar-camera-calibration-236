#pragma once

#include "camcal236/types.hpp"

#include <string>
#include <vector>

namespace camcal236::internal {

// Degeneracy tolerances (see README "Degeneracy tolerances").
// A point set is treated as collinear when the smallest singular value of its
// Hartley-normalised centred coordinates falls below this absolute bound.
constexpr double kCollinearityTol = 1e-6;
// A homography / intrinsic system is treated as singular when the relevant
// smallest singular value or pivot magnitude falls below this bound.
constexpr double kSingularTol = 1e-12;

// Returns empty string on success, otherwise the rejection reason.
std::string Validate(const std::vector<TargetPoint>& target,
                     const std::vector<View>& views, const ImageSize& image);

// Normalised (Hartley) DLT homography from 2D->2D correspondences.
// Returns false when the geometry is degenerate.
bool EstimateHomography(const std::vector<Eigen::Vector2d>& src,
                        const std::vector<Eigen::Vector2d>& dst,
                        Eigen::Matrix3d* h);

// Zhang-style closed-form initialisation. obj[v][i] <-> pixels[v][i] are the
// 2D-2D correspondences of view v. Returns false (with reason) when the
// shared intrinsics are not identifiable from the given homographies.
bool Initialise(const std::vector<std::vector<Eigen::Vector2d>>& obj,
                const std::vector<std::vector<Eigen::Vector2d>>& pixels,
                Intrinsics* k, std::vector<Pose>* poses, std::string* error);

// Joint Levenberg-Marquardt refinement of shared intrinsics and all poses.
// Only cost-non-increasing updates are accepted; fx/fy stay positive and
// rotations stay on SO(3) via Rodrigues updates.
void Refine(const std::vector<std::vector<Eigen::Vector2d>>& obj,
            const std::vector<std::vector<Eigen::Vector2d>>& pixels,
            const Options& options, Intrinsics* k, std::vector<Pose>* poses,
            bool* converged, int* iterations);

}  // namespace camcal236::internal

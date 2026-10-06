#pragma once

#include "camcal236/types.hpp"

namespace camcal236 {

// Projects a target-plane point (millimetres, z = 0) into pixels using the
// shared pinhole model (zero skew, no distortion).
Eigen::Vector2d Project(const Intrinsics& k, const Pose& pose,
                        const Eigen::Vector2d& target_xy_mm);

// Validates inputs, initialises intrinsics + poses from multi-view
// homographies, then jointly refines all parameters by minimising the sum of
// squared 2D pixel residuals over every observation of every view.
//
// On validation or initialisation failure, ok == false and no parameters are
// delivered. If the joint optimiser reaches the iteration cap, ok == true,
// result.converged == false and the last valid (best-cost) estimate is kept.
// Inputs are never modified.
Output Calibrate(const std::vector<TargetPoint>& target,
                 const std::vector<View>& views, const ImageSize& image,
                 const Options& options = Options{});

}  // namespace camcal236

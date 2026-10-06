#include "camcal236/calibration.hpp"

namespace camcal236 {

Eigen::Vector2d Project(const Intrinsics& k, const Pose& pose,
                        const Eigen::Vector2d& target_xy_mm) {
  const Eigen::Vector3d p(target_xy_mm.x(), target_xy_mm.y(), 0.0);
  const Eigen::Vector3d c = pose.R * p + pose.t;
  return {k.fx * (c.x() / c.z()) + k.cx, k.fy * (c.y() / c.z()) + k.cy};
}

}  // namespace camcal236

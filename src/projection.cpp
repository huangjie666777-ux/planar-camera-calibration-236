#include "camcal236/camcal.hpp"

namespace camcal236 {

Eigen::Vector2d project(double fx, double fy, double cx, double cy,
                        const Eigen::Matrix3d& rotation,
                        const Eigen::Vector3d& translation,
                        const Eigen::Vector2d& point_mm) {
    const Eigen::Vector3d X(point_mm.x(), point_mm.y(), 0.0);
    const Eigen::Vector3d c = rotation * X + translation;
    return {fx * c.x() / c.z() + cx, fy * c.y() / c.z() + cy};
}

}  // namespace camcal236


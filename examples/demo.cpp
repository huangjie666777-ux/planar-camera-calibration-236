// Known-parameter multi-view example: synthesizes observations from a ground
// truth camera, calibrates, and compares.
#include "camcal236/camcal.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <random>

using namespace camcal236;

namespace {

Eigen::Matrix3d rot_from_euler(double rx, double ry, double rz) {
    const Eigen::AngleAxisd ax(rx, Eigen::Vector3d::UnitX());
    const Eigen::AngleAxisd ay(ry, Eigen::Vector3d::UnitY());
    const Eigen::AngleAxisd az(rz, Eigen::Vector3d::UnitZ());
    return (az * ay * ax).toRotationMatrix();
}

}  // namespace

int main() {
    // Ground truth intrinsics.
    const double fx = 820.0, fy = 815.0, cx = 640.0, cy = 360.0;
    const int W = 1280, H = 720;

    // 8x6 planar grid, 25 mm spacing (48 target points, IDs 0..47).
    std::vector<TargetPoint> target;
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 8; ++c)
            target.push_back({r * 8 + c, 25.0 * c, 25.0 * r});

    // Four views with distinct orientations.
    const std::vector<Eigen::Matrix3d> R_true = {
        rot_from_euler(0.15, -0.10, 0.05),
        rot_from_euler(-0.25, 0.20, -0.10),
        rot_from_euler(0.30, 0.15, 0.20),
        rot_from_euler(-0.10, -0.30, 0.15),
    };
    const std::vector<Eigen::Vector3d> t_true = {
        {-60.0, -40.0, 420.0},
        {50.0, -55.0, 380.0},
        {-30.0, 45.0, 450.0},
        {70.0, 30.0, 400.0},
    };

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.15);  // 0.15 px noise

    std::vector<ViewInput> views(R_true.size());
    for (std::size_t v = 0; v < R_true.size(); ++v) {
        views[v].width = W;
        views[v].height = H;
        for (const auto& tp : target) {
            // Simulate missing points: skip ~10% of observations.
            if ((tp.id + static_cast<int>(v)) % 11 == 0) continue;
            const Eigen::Vector2d p =
                project(fx, fy, cx, cy, R_true[v], t_true[v], {tp.x_mm, tp.y_mm});
            views[v].observations.push_back(
                {tp.id, p.x() + noise(rng), p.y() + noise(rng)});
        }
    }

    CalibrationResult res = calibrate(target, views);
    if (!res.success) {
        std::printf("calibration failed: %s\n", res.message.c_str());
        return 1;
    }

    std::printf("status: %s (converged=%d)\n", res.message.c_str(),
                static_cast<int>(res.converged));
    std::printf("intrinsics  est: fx=%.3f fy=%.3f cx=%.3f cy=%.3f\n", res.fx,
                res.fy, res.cx, res.cy);
    std::printf("intrinsics true: fx=%.3f fy=%.3f cx=%.3f cy=%.3f\n", fx, fy,
                cx, cy);
    for (std::size_t v = 0; v < res.views.size(); ++v) {
        const auto& vr = res.views[v];
        std::printf("view %zu: rms=%.4f px  t_est=(%.2f %.2f %.2f) t_true=(%.2f %.2f %.2f)\n",
                    v, vr.rms, vr.translation.x(), vr.translation.y(),
                    vr.translation.z(), t_true[v].x(), t_true[v].y(),
                    t_true[v].z());
    }
    std::printf("overall rms = %.4f px over %zu observations\n", res.rms,
                res.points.size());
    return 0;
}


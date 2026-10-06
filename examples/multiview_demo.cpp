// Synthetic multi-view calibration example with known ground truth.
// A 8x6 planar grid (48 points, 10 mm spacing) is observed by 4 views with
// known intrinsics/poses plus Gaussian pixel noise, then calibrated.
#include <cmath>
#include <cstdio>
#include <random>

#include "camcal236/calibration.hpp"

namespace {

Eigen::Matrix3d RotFromEuler(double rx, double ry, double rz) {
  Eigen::Matrix3d x, y, z;
  x << 1, 0, 0, 0, std::cos(rx), -std::sin(rx), 0, std::sin(rx), std::cos(rx);
  y << std::cos(ry), 0, std::sin(ry), 0, 1, 0, -std::sin(ry), 0, std::cos(ry);
  z << std::cos(rz), -std::sin(rz), 0, std::sin(rz), std::cos(rz), 0, 0, 0, 1;
  return z * y * x;
}

}  // namespace

int main() {
  using namespace camcal236;

  const ImageSize image{1280, 960};
  const Intrinsics truth_k{820.0, 815.0, 640.0, 480.0};

  std::vector<TargetPoint> target;
  int id = 0;
  for (int row = 0; row < 6; ++row)
    for (int col = 0; col < 8; ++col)
      target.push_back({id++, col * 10.0, row * 10.0});

  std::vector<Pose> truth_poses(4);
  truth_poses[0] = {RotFromEuler(0.15, -0.10, 0.05), {-30, -20, 320}};
  truth_poses[1] = {RotFromEuler(-0.25, 0.20, -0.10), {-35, -25, 300}};
  truth_poses[2] = {RotFromEuler(0.30, 0.15, 0.20), {-25, -30, 340}};
  truth_poses[3] = {RotFromEuler(-0.10, -0.30, 0.15), {-40, -15, 310}};

  std::mt19937 rng(42);
  std::normal_distribution<double> noise(0.0, 0.3);  // 0.3 px sigma

  std::vector<View> views(truth_poses.size());
  for (size_t v = 0; v < truth_poses.size(); ++v) {
    for (const auto& tp : target) {
      // Simulate a missing point in view 1 to exercise id association.
      if (v == 1 && tp.id == 17) continue;
      const Eigen::Vector2d uv =
          Project(truth_k, truth_poses[v], {tp.x_mm, tp.y_mm});
      views[v].observations.push_back(
          {tp.id, uv.x() + noise(rng), uv.y() + noise(rng)});
    }
  }

  const Output out = Calibrate(target, views, image);
  if (!out.ok) {
    std::fprintf(stderr, "calibration failed: %s\n", out.error.c_str());
    return 1;
  }
  const CalibrationResult& r = out.result;
  std::printf("converged=%d iterations=%d total_rms=%.4f px\n", r.converged,
              r.iterations, r.total_rms);
  std::printf("intrinsics: fx=%.3f (%.1f) fy=%.3f (%.1f) cx=%.3f (%.1f) "
              "cy=%.3f (%.1f)\n",
              r.intrinsics.fx, truth_k.fx, r.intrinsics.fy, truth_k.fy,
              r.intrinsics.cx, truth_k.cx, r.intrinsics.cy, truth_k.cy);
  for (size_t v = 0; v < r.views.size(); ++v) {
    const auto& vr = r.views[v];
    std::printf("view %zu: rms=%.4f px  t_err=%.4f mm  R_err=%.5f\n", v,
                vr.rms, (vr.pose.t - truth_poses[v].t).norm(),
                (vr.pose.R - truth_poses[v].R).norm());
  }
  const double k_err =
      std::abs(r.intrinsics.fx - truth_k.fx) / truth_k.fx +
      std::abs(r.intrinsics.fy - truth_k.fy) / truth_k.fy;
  if (!r.converged || r.total_rms > 0.5 || k_err > 0.02) {
    std::fprintf(stderr, "demo accuracy check failed\n");
    return 1;
  }
  std::printf("demo OK\n");
  return 0;
}

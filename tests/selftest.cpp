// Self-tests: validation rejections, projection, calibration accuracy,
// non-convergence flagging, and input immutability.
#include <cmath>
#include <cstdio>
#include <cstring>

#include "camcal236/calibration.hpp"

namespace {

int failures = 0;

void Check(bool cond, const char* name) {
  std::printf("[%s] %s\n", cond ? "PASS" : "FAIL", name);
  if (!cond) ++failures;
}

Eigen::Matrix3d RotEuler(double rx, double ry, double rz) {
  Eigen::Matrix3d x, y, z;
  x << 1, 0, 0, 0, std::cos(rx), -std::sin(rx), 0, std::sin(rx), std::cos(rx);
  y << std::cos(ry), 0, std::sin(ry), 0, 1, 0, -std::sin(ry), 0, std::cos(ry);
  z << std::cos(rz), -std::sin(rz), 0, std::sin(rz), std::cos(rz), 0, 0, 0, 1;
  return z * y * x;
}

using namespace camcal236;

std::vector<TargetPoint> MakeTarget() {
  std::vector<TargetPoint> t;
  int id = 0;
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 8; ++c) t.push_back({id++, c * 12.0, r * 12.0});
  return t;
}

std::vector<View> MakeViews(const std::vector<TargetPoint>& target,
                            const Intrinsics& k, std::vector<Pose>* poses) {
  *poses = {{RotEuler(0.2, -0.1, 0.05), {-40, -15, 400}},
            {RotEuler(-0.25, 0.2, -0.1), {-45, -20, 380}},
            {RotEuler(0.3, 0.15, 0.2), {-35, -25, 420}}};
  std::vector<View> views(poses->size());
  for (size_t v = 0; v < poses->size(); ++v)
    for (const auto& tp : target) {
      if (v == 2 && tp.id == 5) continue;  // exercise missing points
      const Eigen::Vector2d uv =
          Project(k, (*poses)[v], {tp.x_mm, tp.y_mm});
      views[v].observations.push_back({tp.id, uv.x(), uv.y()});
    }
  return views;
}

}  // namespace

int main() {
  const ImageSize image{1280, 960};
  const Intrinsics k{800.0, 790.0, 640.0, 480.0};
  const auto target = MakeTarget();
  std::vector<Pose> poses;
  const auto views = MakeViews(target, k, &poses);

  // Projection interface: known pinhole geometry.
  {
    Pose p;
    p.t = {0, 0, 100};
    const Eigen::Vector2d uv = Project(k, p, {10, 20});
    Check(std::abs(uv.x() - (640 + 80)) < 1e-9 &&
              std::abs(uv.y() - (480 + 158)) < 1e-9,
          "project: pinhole model");
  }

  // Noiseless calibration recovers ground truth.
  {
    const Output out = Calibrate(target, views, image);
    Check(out.ok, "calibrate: succeeds on clean data");
    Check(out.ok && out.result.converged, "calibrate: converged");
    Check(out.ok && out.result.total_rms < 1e-6, "calibrate: zero residual");
    Check(out.ok && std::abs(out.result.intrinsics.fx - k.fx) < 1e-3 &&
              std::abs(out.result.intrinsics.cx - k.cx) < 1e-3,
          "calibrate: intrinsics recovered");
    Check(out.ok &&
              (out.result.views[0].pose.t - poses[0].t).norm() < 1e-3,
          "calibrate: pose recovered");
    Check(out.ok && out.result.views.size() == 3 &&
              out.result.views[2].points.size() == 31,
          "calibrate: per-point results with missing id");
  }

  // Iteration cap: keeps estimate, flags not converged.
  {
    Options opt;
    opt.max_iterations = 1;
    const Output out = Calibrate(target, views, image, opt);
    Check(out.ok && !out.result.converged && out.result.iterations == 1,
          "calibrate: iteration cap flags non-convergence");
  }

  // Validation rejections.
  {
    auto bad = views;
    bad[0].observations[0].id = 999;
    Check(!Calibrate(target, bad, image).ok, "validate: unknown id");
  }
  {
    auto bad = views;
    bad[0].observations[1].id = bad[0].observations[0].id;
    Check(!Calibrate(target, bad, image).ok, "validate: duplicate observation");
  }
  {
    auto bad = views;
    bad[0].observations[0].u = 1280.0;
    Check(!Calibrate(target, bad, image).ok, "validate: out of image");
  }
  {
    auto bad = views;
    bad[0].observations.resize(7);
    Check(!Calibrate(target, bad, image).ok, "validate: fewer than 8 points");
  }
  {
    Check(!Calibrate(target, {views[0], views[1]}, image).ok,
          "validate: fewer than 3 views");
  }
  {
    auto col = target;
    for (auto& tp : col) tp.y_mm = 2.0 * tp.x_mm;
    Check(!Calibrate(col, views, image).ok, "validate: collinear target");
  }
  {
    auto dup = target;
    dup[1].id = dup[0].id;
    Check(!Calibrate(dup, views, image).ok, "validate: duplicate target id");
  }
  {
    Check(!Calibrate(target, views, {0, 960}).ok,
          "validate: non-positive image size");
  }

  // Inputs are not modified.
  {
    auto copy_target = target;
    auto copy_views = views;
    const Output out = Calibrate(copy_target, copy_views, image);
    bool same = out.ok;
    for (size_t i = 0; i < target.size(); ++i)
      same &= copy_target[i].x_mm == target[i].x_mm &&
              copy_target[i].y_mm == target[i].y_mm;
    for (size_t v = 0; v < views.size(); ++v)
      for (size_t i = 0; i < views[v].observations.size(); ++i)
        same &= copy_views[v].observations[i].u == views[v].observations[i].u;
    Check(same, "inputs unchanged");
  }

  std::printf(failures == 0 ? "ALL TESTS PASSED\n" : "%d FAILURES\n",
              failures);
  return failures == 0 ? 0 : 1;
}

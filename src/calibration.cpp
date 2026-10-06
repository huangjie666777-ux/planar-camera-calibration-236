#include "camcal236/calibration.hpp"

#include <cmath>
#include <map>

#include "internal.hpp"

namespace camcal236 {

Output Calibrate(const std::vector<TargetPoint>& target,
                 const std::vector<View>& views, const ImageSize& image,
                 const Options& options) {
  Output out;
  const std::string err = internal::Validate(target, views, image);
  if (!err.empty()) {
    out.error = err;
    return out;
  }

  std::map<int, Eigen::Vector2d> by_id;
  for (const auto& tp : target)
    by_id.emplace(tp.id, Eigen::Vector2d(tp.x_mm, tp.y_mm));

  // Per-view 2D-2D correspondences (missing ids simply absent per view).
  std::vector<std::vector<Eigen::Vector2d>> pix(views.size());
  std::vector<std::vector<int>> pix_ids(views.size());
  std::vector<std::vector<Eigen::Vector2d>> obj_per_view(views.size());
  for (size_t v = 0; v < views.size(); ++v) {
    for (const auto& ob : views[v].observations) {
      pix[v].emplace_back(ob.u, ob.v);
      pix_ids[v].push_back(ob.id);
      obj_per_view[v].push_back(by_id.at(ob.id));
    }
  }

  Intrinsics k;
  std::vector<Pose> poses;
  std::string init_err;
  if (!internal::Initialise(obj_per_view, pix, &k, &poses, &init_err)) {
    out.error = init_err;
    return out;
  }

  // Joint refinement over shared intrinsics and all view poses.
  bool converged = false;
  int iterations = 0;
  internal::Refine(obj_per_view, pix, options, &k, &poses, &converged,
                   &iterations);

  // Recompute predictions, residuals and RMS from the final parameters.
  out.result.converged = converged;
  out.result.iterations = iterations;
  out.result.intrinsics = k;
  double sum_sq = 0.0;
  size_t count = 0;
  for (size_t v = 0; v < views.size(); ++v) {
    ViewResult vr;
    vr.pose = poses[v];
    double vsum = 0.0;
    for (size_t i = 0; i < pix_ids[v].size(); ++i) {
      PointResult pr;
      pr.id = pix_ids[v][i];
      pr.predicted = Project(k, poses[v], obj_per_view[v][i]);
      pr.residual = pix[v][i] - pr.predicted;
      vsum += pr.residual.squaredNorm();
      vr.points.push_back(pr);
    }
    vr.rms = std::sqrt(vsum / static_cast<double>(vr.points.size()));
    sum_sq += vsum;
    count += vr.points.size();
    out.result.views.push_back(std::move(vr));
  }
  out.result.total_rms = std::sqrt(sum_sq / static_cast<double>(count));
  out.ok = true;
  return out;
}

}  // namespace camcal236

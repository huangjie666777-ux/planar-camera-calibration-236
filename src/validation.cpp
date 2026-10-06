#include "internal.hpp"

#include <cmath>
#include <map>
#include <set>

namespace camcal236::internal {
namespace {

constexpr int kMaxTargetPoints = 48;
constexpr int kMinViews = 3;
constexpr int kMaxViews = 6;
constexpr int kMinObservationsPerView = 8;

bool Finite(double v) { return std::isfinite(v); }

// True when the set is (numerically) collinear or degenerate.
bool Collinear(const std::vector<Eigen::Vector2d>& pts) {
  if (pts.size() < 3) return true;
  Eigen::Vector2d mean = Eigen::Vector2d::Zero();
  for (const auto& p : pts) mean += p;
  mean /= static_cast<double>(pts.size());
  Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
  for (const auto& p : pts) {
    const Eigen::Vector2d d = p - mean;
    cov += d * d.transpose();
  }
  Eigen::JacobiSVD<Eigen::Matrix2d> svd(cov);
  const auto& sv = svd.singularValues();
  const double scale = std::max({sv(0), 1e-300});
  return sv(1) <= kCollinearityTol * kCollinearityTol * scale *
                      static_cast<double>(pts.size());
}

}  // namespace

std::string Validate(const std::vector<TargetPoint>& target,
                     const std::vector<View>& views, const ImageSize& image) {
  if (target.empty() || target.size() > kMaxTargetPoints)
    return "target must contain 1..48 points";
  if (views.size() < kMinViews || views.size() > kMaxViews)
    return "need 3..6 views";
  if (image.width <= 0 || image.height <= 0)
    return "image width/height must be positive integers";

  std::map<int, Eigen::Vector2d> by_id;
  std::vector<Eigen::Vector2d> target_xy;
  for (const auto& tp : target) {
    if (!Finite(tp.x_mm) || !Finite(tp.y_mm))
      return "target point coordinates must be finite";
    if (!by_id.emplace(tp.id, Eigen::Vector2d(tp.x_mm, tp.y_mm)).second)
      return "duplicate target point id " + std::to_string(tp.id);
    target_xy.emplace_back(tp.x_mm, tp.y_mm);
  }
  if (Collinear(target_xy)) return "target points are collinear";

  for (size_t vi = 0; vi < views.size(); ++vi) {
    const View& view = views[vi];
    const std::string tag = "view " + std::to_string(vi) + ": ";
    if (view.observations.size() < kMinObservationsPerView)
      return tag + "needs at least 8 observations";
    std::set<int> seen;
    std::vector<Eigen::Vector2d> pix;
    std::vector<Eigen::Vector2d> obj;
    for (const auto& ob : view.observations) {
      if (!Finite(ob.u) || !Finite(ob.v))
        return tag + "observation coordinates must be finite";
      if (by_id.count(ob.id) == 0)
        return tag + "unknown target id " + std::to_string(ob.id);
      if (!seen.insert(ob.id).second)
        return tag + "duplicate observation of id " + std::to_string(ob.id);
      if (ob.u < 0.0 || ob.u >= image.width || ob.v < 0.0 ||
          ob.v >= image.height)
        return tag + "observation outside image bounds";
      pix.emplace_back(ob.u, ob.v);
      obj.push_back(by_id.at(ob.id));
    }
    if (Collinear(pix))
      return tag + "observed pixels are collinear (degenerate geometry)";
    if (Collinear(obj))
      return tag + "observed target points are collinear";
  }
  return "";
}

}  // namespace camcal236::internal

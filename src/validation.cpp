#include "detail.hpp"

#include <Eigen/SVD>

#include <cmath>
#include <unordered_set>

namespace camcal236::detail {

namespace {

bool finite(double v) { return std::isfinite(v); }

// Spread of a point set along its least-spread in-plane direction, normalized
// by the point scale. Returns the ratio; small means (near-)collinear.
double normalized_min_spread(const std::vector<Eigen::Vector2d>& pts) {
    Eigen::Vector2d mean = Eigen::Vector2d::Zero();
    for (const auto& p : pts) mean += p;
    mean /= static_cast<double>(pts.size());
    Eigen::Matrix2Xd centered(2, pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) centered.col(i) = pts[i] - mean;
    Eigen::JacobiSVD<Eigen::Matrix2Xd> svd(centered, Eigen::ComputeThinU);
    const auto& sv = svd.singularValues();
    if (sv(0) <= 0.0) return 0.0;
    return sv(1) / sv(0);  // scale-invariant spread ratio
}

}  // namespace

std::string validate(const std::vector<TargetPoint>& target_points,
                     const std::vector<ViewInput>& views,
                     const SolverOptions& options,
                     Problem& problem) {
    if (target_points.empty()) return "no target points provided";
    if (target_points.size() > 48) return "more than 48 target points";
    if (views.size() < 3 || views.size() > 6) return "need 3 to 6 views";

    std::unordered_set<int> ids;
    for (const auto& tp : target_points) {
        if (!finite(tp.x_mm) || !finite(tp.y_mm))
            return "target point " + std::to_string(tp.id) + " is not finite";
        if (!ids.insert(tp.id).second)
            return "duplicate target point id " + std::to_string(tp.id);
    }

    // Target points must not be collinear (planar spread required).
    std::vector<Eigen::Vector2d> all;
    all.reserve(target_points.size());
    for (const auto& tp : target_points) all.emplace_back(tp.x_mm, tp.y_mm);
    if (all.size() >= 3 &&
        normalized_min_spread(all) < options.collinearity_tolerance)
        return "target points are (near-)collinear";

    problem.points = target_points;
    problem.index_of.clear();
    for (std::size_t i = 0; i < target_points.size(); ++i)
        problem.index_of[target_points[i].id] = i;
    problem.raw_views = views;
    problem.views.assign(views.size(), {});

    for (std::size_t v = 0; v < views.size(); ++v) {
        const auto& view = views[v];
        if (view.width <= 0 || view.height <= 0)
            return "view " + std::to_string(v) + " has non-positive image size";
        if (view.observations.size() < 8)
            return "view " + std::to_string(v) + " has fewer than 8 observations";

        std::unordered_set<int> seen;
        std::vector<Eigen::Vector2d> used;
        for (const auto& ob : view.observations) {
            if (!finite(ob.u) || !finite(ob.v))
                return "view " + std::to_string(v) + " has a non-finite observation";
            auto it = problem.index_of.find(ob.id);
            if (it == problem.index_of.end())
                return "view " + std::to_string(v) + " observes unknown id " +
                       std::to_string(ob.id);
            if (!seen.insert(ob.id).second)
                return "view " + std::to_string(v) + " has duplicate observation of id " +
                       std::to_string(ob.id);
            if (ob.u < 0.0 || ob.u > view.width - 1.0 || ob.v < 0.0 ||
                ob.v > view.height - 1.0)
                return "view " + std::to_string(v) + " observation of id " +
                       std::to_string(ob.id) + " lies outside the image";
            problem.views[v].emplace_back(it->second, Eigen::Vector2d(ob.u, ob.v));
            used.emplace_back(target_points[it->second].x_mm,
                              target_points[it->second].y_mm);
        }
        if (normalized_min_spread(used) < options.collinearity_tolerance)
            return "view " + std::to_string(v) +
                   " observes only (near-)collinear target points";
    }
    return {};
}

}  // namespace camcal236::detail

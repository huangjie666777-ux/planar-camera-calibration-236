#include "camcal236/camcal.hpp"

#include "detail.hpp"

#include <cmath>

namespace camcal236 {

CalibrationResult calibrate(const std::vector<TargetPoint>& target_points,
                            const std::vector<ViewInput>& views,
                            const SolverOptions& options) {
    CalibrationResult result;

    detail::Problem problem;
    const std::string err = detail::validate(target_points, views, options, problem);
    if (!err.empty()) {
        result.message = "validation failed: " + err;
        return result;
    }

    double fx, fy, cx, cy;
    std::vector<Eigen::Matrix3d> rotations;
    std::vector<Eigen::Vector3d> translations;
    std::string init_err;
    if (!detail::initialize(problem, options, fx, fy, cx, cy, rotations,
                            translations, init_err)) {
        // Initialization failure: deliver no parameters.
        result.message = "initialization failed: " + init_err;
        return result;
    }

    bool converged = false;
    int iterations = 0;
    detail::refine(problem, options, fx, fy, cx, cy, rotations, translations,
                   converged, iterations);

    // Recompute predictions, residuals and RMS from the final parameters.
    result.fx = fx;
    result.fy = fy;
    result.cx = cx;
    result.cy = cy;
    result.views.resize(views.size());
    double sum_sq = 0.0;
    std::size_t count = 0;
    for (std::size_t v = 0; v < views.size(); ++v) {
        result.views[v].rotation = rotations[v];
        result.views[v].translation = translations[v];
        double view_sq = 0.0;
        for (const auto& [idx, px] : problem.views[v]) {
            const Eigen::Vector2d pred =
                project(fx, fy, cx, cy, rotations[v], translations[v],
                        {problem.points[idx].x_mm, problem.points[idx].y_mm});
            const Eigen::Vector2d res = px - pred;
            PointResult pr;
            pr.view_index = static_cast<int>(v);
            pr.id = problem.points[idx].id;
            pr.predicted = pred;
            pr.residual = res;
            result.points.push_back(pr);
            view_sq += res.squaredNorm();
        }
        result.views[v].rms =
            std::sqrt(view_sq / static_cast<double>(problem.views[v].size()));
        sum_sq += view_sq;
        count += problem.views[v].size();
    }
    result.rms = std::sqrt(sum_sq / static_cast<double>(count));

    result.success = true;
    result.converged = converged;
    result.message = converged
                         ? "converged"
                         : "iteration limit reached; returning last valid estimate";
    return result;
}

}  // namespace camcal236


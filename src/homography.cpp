#include "detail.hpp"

#include <Eigen/SVD>

namespace camcal236::detail {

namespace {

// Hartley isotropic normalization: centroid to origin, mean distance sqrt(2).
Eigen::Matrix3d normalizer(const std::vector<Eigen::Vector2d>& pts) {
    Eigen::Vector2d mean = Eigen::Vector2d::Zero();
    for (const auto& p : pts) mean += p;
    mean /= static_cast<double>(pts.size());
    double mean_dist = 0.0;
    for (const auto& p : pts) mean_dist += (p - mean).norm();
    mean_dist /= static_cast<double>(pts.size());
    const double s = mean_dist > 0.0 ? std::sqrt(2.0) / mean_dist : 1.0;
    Eigen::Matrix3d T;
    T << s, 0, -s * mean.x(), 0, s, -s * mean.y(), 0, 0, 1;
    return T;
}

}  // namespace

bool estimate_homography(const Problem& problem, std::size_t view,
                         Eigen::Matrix3d& H) {
    const auto& obs = problem.views[view];
    const std::size_t n = obs.size();
    std::vector<Eigen::Vector2d> src, dst;
    src.reserve(n);
    dst.reserve(n);
    for (const auto& [idx, px] : obs) {
        src.emplace_back(problem.points[idx].x_mm, problem.points[idx].y_mm);
        dst.push_back(px);
    }
    const Eigen::Matrix3d Ts = normalizer(src);
    const Eigen::Matrix3d Td = normalizer(dst);

    Eigen::MatrixXd A(2 * n, 9);
    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d X = Ts * src[i].homogeneous();
        const Eigen::Vector3d x = Td * dst[i].homogeneous();
        A.row(2 * i) << 0, 0, 0, -X.x(), -X.y(), -X.z(), x.y() * X.x(),
            x.y() * X.y(), x.y() * X.z();
        A.row(2 * i + 1) << X.x(), X.y(), X.z(), 0, 0, 0, -x.x() * X.x(),
            -x.x() * X.y(), -x.x() * X.z();
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
    if (svd.singularValues()(7) <= 1e-12) return false;  // degenerate geometry
    const Eigen::VectorXd h = svd.matrixV().col(8);
    Eigen::Matrix3d Hn;
    Hn << h(0), h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8);
    H = Td.inverse() * Hn * Ts;
    H /= H(2, 2) != 0.0 ? H(2, 2) : 1.0;
    return std::isfinite(H.sum());
}

}  // namespace camcal236::detail


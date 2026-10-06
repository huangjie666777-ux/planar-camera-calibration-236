#include "detail.hpp"

#include <Eigen/SVD>

#include <cmath>

namespace camcal236::detail {

namespace {

// Zhang's v_ij constraint row from homography columns i, j (0-based).
Eigen::Matrix<double, 1, 6> v_row(const Eigen::Matrix3d& H, int i, int j) {
    const auto h = [&](int r, int c) { return H(r, c); };
    Eigen::Matrix<double, 1, 6> v;
    v << h(0, i) * h(0, j),
        h(0, i) * h(1, j) + h(1, i) * h(0, j),
        h(1, i) * h(1, j),
        h(2, i) * h(0, j) + h(0, i) * h(2, j),
        h(2, i) * h(1, j) + h(1, i) * h(2, j),
        h(2, i) * h(2, j);
    return v;
}

Eigen::Matrix3d project_to_so3(const Eigen::Matrix3d& M) {
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
    if (R.determinant() < 0.0) {
        Eigen::Matrix3d U = svd.matrixU();
        U.col(2) *= -1.0;
        R = U * svd.matrixV().transpose();
    }
    return R;
}

}  // namespace

bool initialize(const Problem& problem, const SolverOptions& options,
                double& fx, double& fy, double& cx, double& cy,
                std::vector<Eigen::Matrix3d>& rotations,
                std::vector<Eigen::Vector3d>& translations,
                std::string& error) {
    const std::size_t nv = problem.views.size();
    std::vector<Eigen::Matrix3d> H(nv);
    for (std::size_t v = 0; v < nv; ++v) {
        if (!estimate_homography(problem, v, H[v])) {
            error = "degenerate homography geometry in view " + std::to_string(v);
            return false;
        }
    }

    // Intrinsics from orthogonalities of B = K^-T K^-1 (skew = 0).
    Eigen::MatrixXd V(2 * nv, 6);
    for (std::size_t v = 0; v < nv; ++v) {
        V.row(2 * v) = v_row(H[v], 0, 1);
        V.row(2 * v + 1) = v_row(H[v], 0, 0) - v_row(H[v], 1, 1);
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeFullV);
    const auto& sv = svd.singularValues();
    // b must lie in a 1-dimensional null space of V. If the second-smallest
    // singular value also vanishes, B is not uniquely determined and the
    // intrinsics are unidentifiable (e.g. near-parallel target planes).
    if (sv(4) <= 1e-10 * sv(0)) {
        error = "intrinsics not identifiable: views are too similar "
                "(near-parallel target planes)";
        return false;
    }
    const Eigen::VectorXd b = svd.matrixV().col(5);
    const double B00 = b(0), B01 = b(1), B11 = b(2), B02 = b(3), B12 = b(4),
                 B22 = b(5);

    if (std::abs(B00 * B11 - B01 * B01) < 1e-30) {
        error = "intrinsics not identifiable (singular B matrix)";
        return false;
    }
    const double sign = (B00 > 0.0 && (B00 * B22 - B02 * B02) > 0.0) ? 1.0 : -1.0;
    const double sB00 = sign * B00, sB01 = sign * B01, sB11 = sign * B11,
                 sB02 = sign * B02, sB12 = sign * B12, sB22 = sign * B22;
    const double det = sB00 * sB11 - sB01 * sB01;
    if (sB00 <= 0.0 || det <= 0.0) {
        error = "intrinsics not identifiable: no positive focal solution";
        return false;
    }
    cy = (sB01 * sB02 - sB00 * sB12) / det;
    cx = -sB02 / sB00;
    const double lambda = sB22 - (sB02 * sB02 + cy * cy * det) / sB00;
    if (lambda <= 0.0) {
        error = "intrinsics not identifiable: inconsistent scale";
        return false;
    }
    const double fx2 = lambda / sB00;
    const double fy2 = lambda * sB00 / det;
    if (fx2 <= 0.0 || fy2 <= 0.0) {
        error = "intrinsics not identifiable: non-positive focal length";
        return false;
    }
    fx = std::sqrt(fx2);
    fy = std::sqrt(fy2);

    Eigen::Matrix3d K;
    K << fx, 0, cx, 0, fy, cy, 0, 0, 1;
    const Eigen::Matrix3d Kinv = K.inverse();

    rotations.assign(nv, Eigen::Matrix3d::Identity());
    translations.assign(nv, Eigen::Vector3d::Zero());
    for (std::size_t v = 0; v < nv; ++v) {
        const Eigen::Vector3d h1 = H[v].col(0), h2 = H[v].col(1), h3 = H[v].col(2);
        const double s1 = (Kinv * h1).norm();
        const double s2 = (Kinv * h2).norm();
        if (s1 <= 0.0 || s2 <= 0.0) {
            error = "degenerate homography scale in view " + std::to_string(v);
            return false;
        }
        const double scale = 2.0 / (s1 + s2);
        Eigen::Vector3d r1 = scale * (Kinv * h1);
        Eigen::Vector3d r2 = scale * (Kinv * h2);
        Eigen::Vector3d r3 = r1.cross(r2);
        Eigen::Matrix3d R;
        R.col(0) = r1;
        R.col(1) = r2;
        R.col(2) = r3;
        R = project_to_so3(R);
        Eigen::Vector3d t = scale * (Kinv * h3);

        // Pick the sign that puts the observed target points in front.
        auto mean_depth = [&](const Eigen::Matrix3d& Rr, const Eigen::Vector3d& tt) {
            double z = 0.0;
            for (const auto& [idx, px] : problem.views[v]) {
                const Eigen::Vector3d X(problem.points[idx].x_mm,
                                        problem.points[idx].y_mm, 0.0);
                z += (Rr * X + tt).z();
            }
            return z / static_cast<double>(problem.views[v].size());
        };
        if (mean_depth(R, t) < options.depth_epsilon) {
            R.col(0) = -R.col(0);
            R.col(1) = -R.col(1);
            t = -t;
            if (mean_depth(R, t) < options.depth_epsilon) {
                error = "no positive-depth solution in view " + std::to_string(v);
                return false;
            }
        }
        rotations[v] = R;
        translations[v] = t;
    }
    return true;
}

}  // namespace camcal236::detail

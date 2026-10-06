#include "detail.hpp"

#include <cmath>

namespace camcal236::detail {

namespace {

Eigen::Matrix3d exp_so3(const Eigen::Vector3d& w) {
    const double th = w.norm();
    Eigen::Matrix3d W;
    W << 0, -w.z(), w.y(), w.z(), 0, -w.x(), -w.y(), w.x(), 0;
    if (th < 1e-14) return Eigen::Matrix3d::Identity() + W;
    const double a = std::sin(th) / th;
    const double b = (1.0 - std::cos(th)) / (th * th);
    return Eigen::Matrix3d::Identity() + a * W + b * W * W;
}

struct State {
    double fx, fy, cx, cy;
    std::vector<Eigen::Matrix3d> R;
    std::vector<Eigen::Vector3d> t;
};

// Sum of squared 2D residuals over all observations of all views.
double cost(const Problem& p, const State& s) {
    double c = 0.0;
    for (std::size_t v = 0; v < p.views.size(); ++v) {
        for (const auto& [idx, px] : p.views[v]) {
            const Eigen::Vector3d X(p.points[idx].x_mm, p.points[idx].y_mm, 0.0);
            const Eigen::Vector3d q = s.R[v] * X + s.t[v];
            const double du = s.fx * q.x() / q.z() + s.cx - px.x();
            const double dv = s.fy * q.y() / q.z() + s.cy - px.y();
            c += du * du + dv * dv;
        }
    }
    return c;
}

}  // namespace

void refine(const Problem& problem, const SolverOptions& options,
            double& fx, double& fy, double& cx, double& cy,
            std::vector<Eigen::Matrix3d>& rotations,
            std::vector<Eigen::Vector3d>& translations,
            bool& converged, int& iterations) {
    const std::size_t nv = problem.views.size();
    const std::size_t nparam = 4 + 6 * nv;
    std::size_t nres = 0;
    for (const auto& v : problem.views) nres += 2 * v.size();

    State s{fx, fy, cx, cy, rotations, translations};
    double cur = cost(problem, s);
    double lambda = 1e-3;
    converged = false;
    iterations = 0;

    Eigen::VectorXd r(nres);
    Eigen::MatrixXd J(nres, nparam);

    for (int it = 0; it < options.max_iterations; ++it) {
        iterations = it + 1;
        // Assemble residuals and analytic Jacobian.
        std::size_t row = 0;
        for (std::size_t v = 0; v < nv; ++v) {
            for (const auto& [idx, px] : problem.views[v]) {
                const Eigen::Vector3d X(problem.points[idx].x_mm,
                                        problem.points[idx].y_mm, 0.0);
                const Eigen::Vector3d q = s.R[v] * X + s.t[v];
                const double z = q.z();
                const double zu = s.fx * q.x() / z, zv = s.fy * q.y() / z;
                r(row) = zu + s.cx - px.x();
                r(row + 1) = zv + s.cy - px.y();

                J.row(row).setZero();
                J.row(row + 1).setZero();
                // Intrinsics: fx, fy parameterized in log-space to stay positive.
                J(row, 0) = zu;
                J(row + 1, 1) = zv;
                J(row, 2) = 1.0;
                J(row + 1, 3) = 1.0;
                // Pose of view v: dres/dq * dq/d(omega, t). A left
                // perturbation R <- exp(w) R rotates about the camera
                // center, so dq/dw = -[R X]x (not -[q]x).
                Eigen::Matrix<double, 2, 3> D;
                D << s.fx / z, 0, -s.fx * q.x() / (z * z), 0, s.fy / z,
                    -s.fy * q.y() / (z * z);
                const Eigen::Vector3d rx = s.R[v] * X;
                Eigen::Matrix3d skew;
                skew << 0, -rx.z(), rx.y(), rx.z(), 0, -rx.x(), -rx.y(), rx.x(), 0;
                const Eigen::Matrix<double, 2, 3> dw = -D * skew;  // left perturbation
                const std::size_t base = 4 + 6 * v;
                J.block<1, 3>(row, base) = dw.row(0);
                J.block<1, 3>(row + 1, base) = dw.row(1);
                J.block<1, 3>(row, base + 3) = D.row(0);
                J.block<1, 3>(row + 1, base + 3) = D.row(1);
                row += 2;
            }
        }

        const Eigen::MatrixXd H = J.transpose() * J;
        const Eigen::VectorXd g = J.transpose() * r;
        const double ginf = g.cwiseAbs().maxCoeff();
        if (ginf < 1e-8) {  // gradient vanished: at a stationary point
            converged = true;
            break;
        }

        bool accepted = false;
        double next = cur;
        State trial = s;
        Eigen::VectorXd delta;
        for (int attempt = 0; attempt < 30; ++attempt) {
            Eigen::MatrixXd A = H;
            for (std::size_t k = 0; k < nparam; ++k)
                A(k, k) += lambda * std::max(H(k, k), 1e-12);
            delta = A.ldlt().solve(-g);
            if (!delta.allFinite()) {
                lambda *= 10.0;
                continue;
            }
            trial = s;
            trial.fx = s.fx * std::exp(delta(0));
            trial.fy = s.fy * std::exp(delta(1));
            trial.cx = s.cx + delta(2);
            trial.cy = s.cy + delta(3);
            for (std::size_t v = 0; v < nv; ++v) {
                const Eigen::Vector3d w = delta.segment<3>(4 + 6 * v);
                const Eigen::Vector3d dt = delta.segment<3>(4 + 6 * v + 3);
                trial.R[v] = exp_so3(w) * s.R[v];  // stays on SO(3)
                trial.t[v] = s.t[v] + dt;
            }
            next = cost(problem, trial);
            if (next <= cur) {  // only non-increasing updates are accepted
                accepted = true;
                break;
            }
            lambda *= 10.0;
        }
        if (!accepted) {  // cannot improve; keep last valid estimate
            if (ginf < 1e-6) converged = true;
            break;
        }

        const double rel = (cur - next) / std::max(cur, 1e-30);
        s = trial;
        cur = next;
        lambda = std::max(lambda / 5.0, 1e-12);
        if (rel < options.convergence_tolerance) {
            converged = true;
            break;
        }
    }

    fx = s.fx;
    fy = s.fy;
    cx = s.cx;
    cy = s.cy;
    rotations = s.R;
    translations = s.t;
}

}  // namespace camcal236::detail

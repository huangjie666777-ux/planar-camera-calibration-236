#include "internal.hpp"

#include <cmath>

#include "camcal236/calibration.hpp"

namespace camcal236::internal {
namespace {

Eigen::Vector3d RotationToRodrigues(const Eigen::Matrix3d& r) {
  const double cos_theta = std::clamp((r.trace() - 1.0) * 0.5, -1.0, 1.0);
  const double theta = std::acos(cos_theta);
  if (theta < 1e-12) return Eigen::Vector3d::Zero();
  Eigen::Vector3d axis(r(2, 1) - r(1, 2), r(0, 2) - r(2, 0),
                       r(1, 0) - r(0, 1));
  const double s = 2.0 * std::sin(theta);
  if (std::abs(s) > 1e-12) {
    axis /= s;
  } else {
    // theta near pi: recover axis from the symmetric part of R.
    axis = Eigen::Vector3d(std::sqrt(std::max(0.0, (r(0, 0) + 1.0) * 0.5)),
                           std::sqrt(std::max(0.0, (r(1, 1) + 1.0) * 0.5)),
                           std::sqrt(std::max(0.0, (r(2, 2) + 1.0) * 0.5)));
    if (r(0, 1) < 0) axis.y() = -axis.y();
    if (r(0, 2) < 0) axis.z() = -axis.z();
    if (axis.norm() < 1e-12) axis = Eigen::Vector3d(1, 0, 0);
    axis.normalize();
  }
  return axis * theta;
}

Eigen::Matrix3d RodriguesToRotation(const Eigen::Vector3d& w) {
  const double theta = w.norm();
  Eigen::Matrix3d k;
  k << 0, -w.z(), w.y(), w.z(), 0, -w.x(), -w.y(), w.x(), 0;
  if (theta < 1e-12) return Eigen::Matrix3d::Identity() + k;
  return Eigen::Matrix3d::Identity() + std::sin(theta) / theta * k +
         (1.0 - std::cos(theta)) / (theta * theta) * k * k;
}

struct Problem {
  const std::vector<std::vector<Eigen::Vector2d>>& obj;
  const std::vector<std::vector<Eigen::Vector2d>>& pix;
  int nviews;
  int nparams;
  int nresiduals;
};

void Unpack(const Problem& prob, const Eigen::VectorXd& p, Intrinsics* k,
            std::vector<Pose>* poses) {
  k->fx = p(0);
  k->fy = p(1);
  k->cx = p(2);
  k->cy = p(3);
  poses->resize(prob.nviews);
  for (int v = 0; v < prob.nviews; ++v) {
    (*poses)[v].R = RodriguesToRotation(p.segment<3>(4 + 6 * v));
    (*poses)[v].t = p.segment<3>(7 + 6 * v);
  }
}

Eigen::VectorXd Residuals(const Problem& prob, const Eigen::VectorXd& p) {
  Intrinsics k;
  std::vector<Pose> poses;
  Unpack(prob, p, &k, &poses);
  Eigen::VectorXd r(prob.nresiduals);
  int row = 0;
  for (int v = 0; v < prob.nviews; ++v) {
    for (size_t i = 0; i < prob.pix[v].size(); ++i) {
      const Eigen::Vector2d pred = Project(k, poses[v], prob.obj[v][i]);
      r(row++) = pred.x() - prob.pix[v][i].x();
      r(row++) = pred.y() - prob.pix[v][i].y();
    }
  }
  return r;
}

}  // namespace

void Refine(const std::vector<std::vector<Eigen::Vector2d>>& obj,
            const std::vector<std::vector<Eigen::Vector2d>>& pixels,
            const Options& options, Intrinsics* k, std::vector<Pose>* poses,
            bool* converged, int* iterations) {
  Problem prob{obj, pixels, static_cast<int>(pixels.size()), 0, 0};
  prob.nparams = 4 + 6 * prob.nviews;
  int nres = 0;
  for (const auto& v : pixels) nres += 2 * static_cast<int>(v.size());
  prob.nresiduals = nres;

  Eigen::VectorXd p(prob.nparams);
  p(0) = k->fx;
  p(1) = k->fy;
  p(2) = k->cx;
  p(3) = k->cy;
  for (int v = 0; v < prob.nviews; ++v) {
    p.segment<3>(4 + 6 * v) = RotationToRodrigues((*poses)[v].R);
    p.segment<3>(7 + 6 * v) = (*poses)[v].t;
  }

  Eigen::VectorXd r = Residuals(prob, p);
  double cost = 0.5 * r.squaredNorm();
  double lambda = 1e-3;
  *converged = false;
  *iterations = 0;

  for (int it = 0; it < options.max_iterations; ++it) {
    *iterations = it + 1;

    // Central-difference Jacobian.
    Eigen::MatrixXd j(prob.nresiduals, prob.nparams);
    for (int c = 0; c < prob.nparams; ++c) {
      const double h = 1e-7 * std::max(1.0, std::abs(p(c)));
      Eigen::VectorXd pp = p, pm = p;
      pp(c) += h;
      pm(c) -= h;
      j.col(c) = (Residuals(prob, pp) - Residuals(prob, pm)) / (2.0 * h);
    }
    const Eigen::VectorXd g = j.transpose() * r;
    if (g.cwiseAbs().maxCoeff() < options.convergence_tol) {
      *converged = true;
      break;
    }

    const Eigen::MatrixXd jtj = j.transpose() * j;
    bool stepped = false;
    for (int attempt = 0; attempt < 50; ++attempt) {
      Eigen::MatrixXd a = jtj;
      a.diagonal() += lambda * jtj.diagonal().cwiseMax(1e-12);
      const Eigen::VectorXd delta = a.ldlt().solve(-g);
      if (!delta.allFinite()) {
        lambda *= 10.0;
        continue;
      }
      const Eigen::VectorXd trial = p + delta;
      if (trial(0) <= 0.0 || trial(1) <= 0.0) {  // fx, fy must stay positive
        lambda *= 10.0;
        continue;
      }
      const Eigen::VectorXd rt = Residuals(prob, trial);
      const double trial_cost = 0.5 * rt.squaredNorm();
      if (trial_cost < cost) {  // only accept non-increasing updates
        const double rel = (cost - trial_cost) / std::max(cost, 1e-300);
        const double step = delta.norm() / std::max(p.norm(), 1e-300);
        p = trial;
        r = rt;
        cost = trial_cost;
        lambda = std::max(lambda * 0.3, 1e-12);
        stepped = true;
        if (rel < options.convergence_tol && step < std::sqrt(options.convergence_tol)) {
          *converged = true;
        }
        break;
      }
      lambda *= 10.0;
    }
    if (*converged || !stepped) {
      if (!stepped) *converged = true;  // no improving step exists: stationary
      break;
    }
  }

  Unpack(prob, p, k, poses);
}

}  // namespace camcal236::internal

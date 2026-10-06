#include "internal.hpp"

#include <cmath>

namespace camcal236::internal {
namespace {

Eigen::Matrix3d ProjectToSO3(const Eigen::Matrix3d& m) {
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(m, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix3d r = svd.matrixU() * svd.matrixV().transpose();
  if (r.determinant() < 0) {
    Eigen::Matrix3d u = svd.matrixU();
    u.col(2) *= -1.0;
    r = u * svd.matrixV().transpose();
  }
  return r;
}

// Row of the Zhang V matrix for homography columns i, j (0-based).
Eigen::Matrix<double, 1, 6> VRow(const Eigen::Matrix3d& h, int i, int j) {
  Eigen::Matrix<double, 1, 6> v;
  v << h(0, i) * h(0, j), h(0, i) * h(1, j) + h(1, i) * h(0, j),
      h(1, i) * h(1, j), h(2, i) * h(0, j) + h(0, i) * h(2, j),
      h(2, i) * h(1, j) + h(1, i) * h(2, j), h(2, i) * h(2, j);
  return v;
}

}  // namespace

bool Initialise(const std::vector<std::vector<Eigen::Vector2d>>& obj,
                const std::vector<std::vector<Eigen::Vector2d>>& pixels,
                Intrinsics* k, std::vector<Pose>* poses, std::string* error) {
  const int nviews = static_cast<int>(pixels.size());

  std::vector<Eigen::Matrix3d> homographies(nviews);
  for (int v = 0; v < nviews; ++v) {
    if (!EstimateHomography(obj[v], pixels[v], &homographies[v])) {
      *error = "view " + std::to_string(v) +
               ": homography estimation failed (degenerate geometry)";
      return false;
    }
  }

  // Orthogonality constraints on B = K^-T K^-1 from all views jointly.
  Eigen::MatrixXd vmat(2 * nviews, 6);
  for (int v = 0; v < nviews; ++v) {
    vmat.row(2 * v) = VRow(homographies[v], 0, 1);
    vmat.row(2 * v + 1) =
        VRow(homographies[v], 0, 0) - VRow(homographies[v], 1, 1);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(vmat, Eigen::ComputeFullV);
  // Rank must be exactly 5: a vanishing second-smallest singular value means
  // the null space (hence the intrinsics) is not uniquely determined.
  if (svd.singularValues()(4) < kSingularTol) {
    *error = "intrinsics not identifiable: constraint system is singular "
             "(degenerate view configuration)";
    return false;
  }
  const Eigen::Matrix<double, 6, 1> b = svd.matrixV().col(5);
  const double b11 = b(0), b12 = b(1), b22 = b(2), b13 = b(3), b23 = b(4),
               b33 = b(5);

  // B = K^-T K^-1 must be positive definite; its entries scale with 1/f^2 so
  // only sign (not absolute magnitude) is a valid degeneracy test here.
  const double det = b11 * b22 - b12 * b12;
  if (b11 <= 0.0 || det <= 0.0) {
    *error = "intrinsics not identifiable: singular image-of-absolute-conic";
    return false;
  }
  const double cy = (b12 * b13 - b11 * b23) / det;
  const double lambda = b33 - (b13 * b13 + cy * (b12 * b13 - b11 * b23)) / b11;
  if (lambda <= 0.0 || b11 <= 0.0 || det * b11 <= 0.0) {
    *error = "intrinsics not identifiable: non-positive focal length "
             "(degenerate view configuration)";
    return false;
  }
  const double fx = std::sqrt(lambda / b11);
  const double fy = std::sqrt(lambda * b11 / det);
  const double cx = -(b13 * fx * fx) / lambda;  // skew term fixed to zero
  if (!std::isfinite(fx) || !std::isfinite(fy) || !std::isfinite(cx) ||
      !std::isfinite(cy) || fx <= 0.0 || fy <= 0.0) {
    *error = "intrinsics not identifiable: non-finite result";
    return false;
  }
  k->fx = fx;
  k->fy = fy;
  k->cx = cx;
  k->cy = cy;

  Eigen::Matrix3d km;
  km << fx, 0, cx, 0, fy, cy, 0, 0, 1;
  const Eigen::Matrix3d k_inv = km.inverse();

  poses->assign(nviews, Pose{});
  for (int v = 0; v < nviews; ++v) {
    const Eigen::Matrix3d& h = homographies[v];
    const Eigen::Vector3d h1 = h.col(0), h2 = h.col(1), h3 = h.col(2);
    const double n1 = (k_inv * h1).norm(), n2 = (k_inv * h2).norm();
    if (n1 < kSingularTol || n2 < kSingularTol) {
      *error = "view " + std::to_string(v) + ": degenerate homography scale";
      return false;
    }
    const double scale = 2.0 / (n1 + n2);
    Eigen::Vector3d r1 = scale * (k_inv * h1);
    Eigen::Vector3d r2 = scale * (k_inv * h2);
    Eigen::Vector3d t = scale * (k_inv * h3);

    // Disambiguate the twofold sign ambiguity: keep the solution where the
    // target points lie in front of the camera (positive depth).
    double zsum = 0.0;
    for (const auto& p : obj[v])
      zsum += r1.z() * p.x() + r2.z() * p.y() + t.z();
    if (zsum < 0.0) {
      r1 = -r1;
      r2 = -r2;
      t = -t;
    }
    Eigen::Matrix3d r;
    r.col(0) = r1;
    r.col(1) = r2;
    r.col(2) = r1.cross(r2);
    (*poses)[v].R = ProjectToSO3(r);
    (*poses)[v].t = t;
  }
  return true;
}

}  // namespace camcal236::internal

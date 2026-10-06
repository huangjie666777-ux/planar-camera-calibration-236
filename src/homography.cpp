#include "internal.hpp"

namespace camcal236::internal {
namespace {

// Hartley normalisation: translate to centroid, scale mean distance to sqrt(2).
Eigen::Matrix3d Normaliser(const std::vector<Eigen::Vector2d>& pts,
                           std::vector<Eigen::Vector2d>* out) {
  Eigen::Vector2d mean = Eigen::Vector2d::Zero();
  for (const auto& p : pts) mean += p;
  mean /= static_cast<double>(pts.size());
  double mean_dist = 0.0;
  for (const auto& p : pts) mean_dist += (p - mean).norm();
  mean_dist /= static_cast<double>(pts.size());
  const double s =
      (mean_dist > kSingularTol) ? std::sqrt(2.0) / mean_dist : 1.0;
  out->reserve(pts.size());
  for (const auto& p : pts) out->push_back(s * (p - mean));
  Eigen::Matrix3d t;
  t << s, 0, -s * mean.x(), 0, s, -s * mean.y(), 0, 0, 1;
  return t;
}

}  // namespace

bool EstimateHomography(const std::vector<Eigen::Vector2d>& src,
                        const std::vector<Eigen::Vector2d>& dst,
                        Eigen::Matrix3d* h) {
  if (src.size() != dst.size() || src.size() < 4) return false;
  std::vector<Eigen::Vector2d> sn, dn;
  const Eigen::Matrix3d ts = Normaliser(src, &sn);
  const Eigen::Matrix3d td = Normaliser(dst, &dn);

  const int n = static_cast<int>(src.size());
  Eigen::MatrixXd a(2 * n, 9);
  for (int i = 0; i < n; ++i) {
    const double x = sn[i].x(), y = sn[i].y();
    const double u = dn[i].x(), v = dn[i].y();
    a.row(2 * i) << -x, -y, -1, 0, 0, 0, u * x, u * y, u;
    a.row(2 * i + 1) << 0, 0, 0, -x, -y, -1, v * x, v * y, v;
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(a, Eigen::ComputeFullV);
  // Rank must be exactly 8: a vanishing second-smallest singular value means
  // the correspondences do not constrain a unique homography.
  if (svd.singularValues()(7) < kSingularTol) return false;
  Eigen::Matrix3d hn;
  hn.row(0) = svd.matrixV().col(8).segment<3>(0);
  hn.row(1) = svd.matrixV().col(8).segment<3>(3);
  hn.row(2) = svd.matrixV().col(8).segment<3>(6);
  *h = td.inverse() * hn * ts;
  if (std::abs((*h)(2, 2)) > kSingularTol) *h /= (*h)(2, 2);
  return h->allFinite();
}

}  // namespace camcal236::internal

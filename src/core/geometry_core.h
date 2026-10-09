// core/geometry_core.h
// Planar geometry for area along an accumulation: convex hull area of the
// first k sites of an order (incremental), and Voronoi cell areas of all
// sites clipped to a rectangle. Pure C++, coordinates through x(i) / y(i).

#ifndef SPACC_CORE_GEOMETRY_H
#define SPACC_CORE_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace spacc {

struct Point2 { double x, y; };

inline double cross(const Point2& o, const Point2& a, const Point2& b) {
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Shoelace area of a simple polygon given in order.
inline double polygon_area(const std::vector<Point2>& poly) {
  double a = 0.0;
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; i++) {
    const Point2& p = poly[i];
    const Point2& q = poly[(i + 1) % n];
    a += p.x * q.y - q.x * p.y;
  }
  return std::abs(a) / 2.0;
}

// Convex hull (Andrew's monotone chain), counter-clockwise, without
// collinear points.
inline std::vector<Point2> convex_hull(std::vector<Point2> pts) {
  std::sort(pts.begin(), pts.end(), [](const Point2& a, const Point2& b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  });
  pts.erase(std::unique(pts.begin(), pts.end(), [](const Point2& a, const Point2& b) {
    return a.x == b.x && a.y == b.y;
  }), pts.end());
  if (pts.size() < 3) return pts;
  std::vector<Point2> hull(2 * pts.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < pts.size(); i++) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0) k--;
    hull[k++] = pts[i];
  }
  for (std::size_t i = pts.size() - 1, t = k + 1; i > 0; i--) {
    while (k >= t && cross(hull[k - 2], hull[k - 1], pts[i - 1]) <= 0) k--;
    hull[k++] = pts[i - 1];
  }
  hull.resize(k - 1);
  return hull;
}

// Point inside or on a counter-clockwise convex polygon with >= 3 vertices.
inline bool in_convex(const std::vector<Point2>& hull, const Point2& p) {
  for (std::size_t i = 0; i < hull.size(); i++) {
    if (cross(hull[i], hull[(i + 1) % hull.size()], p) < 0) return false;
  }
  return true;
}

// Convex hull area of the first k sites along an order, for k = 1..n_steps.
// A site inside the current hull leaves it unchanged; otherwise the hull is
// rebuilt from its vertices and the new site.
template<class XFn, class YFn, class OrderFn>
void hull_area_along_order(const XFn& x, const YFn& y, const OrderFn& order, int n_steps,
                           double* area) {
  std::vector<Point2> hull;
  for (int k = 0; k < n_steps; k++) {
    const int site = order(k);
    const Point2 p{x(site), y(site)};
    if (hull.size() < 3 || !in_convex(hull, p)) {
      hull.push_back(p);
      hull = convex_hull(hull);
    }
    area[k] = hull.size() >= 3 ? polygon_area(hull) : 0.0;
  }
}

// Clip a convex polygon to the half-plane of points closer to a than to b
// (Sutherland-Hodgman against the perpendicular bisector of a and b).
inline std::vector<Point2> clip_to_bisector(const std::vector<Point2>& poly,
                                            const Point2& a, const Point2& b) {
  const double nx = b.x - a.x, ny = b.y - a.y;
  const double c = (b.x * b.x + b.y * b.y - a.x * a.x - a.y * a.y) / 2.0;
  auto side = [&](const Point2& p) { return nx * p.x + ny * p.y - c; };   // <= 0 keeps
  std::vector<Point2> out;
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; i++) {
    const Point2& p = poly[i];
    const Point2& q = poly[(i + 1) % n];
    const double sp = side(p), sq = side(q);
    if (sp <= 0) out.push_back(p);
    if ((sp < 0 && sq > 0) || (sp > 0 && sq < 0)) {
      const double t = sp / (sp - sq);
      out.push_back(Point2{p.x + t * (q.x - p.x), p.y + t * (q.y - p.y)});
    }
  }
  return out;
}

// Area of each site's Voronoi cell within the rectangle [x0, x1] x [y0, y1].
// Sites at identical coordinates share their cell equally.
template<class XFn, class YFn>
std::vector<double> voronoi_cell_areas(const XFn& x, const YFn& y, int n,
                                       double x0, double x1, double y0, double y1) {
  std::vector<double> area(n, 0.0);
  for (int i = 0; i < n; i++) {
    const Point2 a{x(i), y(i)};
    std::vector<Point2> cell = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    int duplicates = 1;
    for (int j = 0; j < n && !cell.empty(); j++) {
      if (j == i) continue;
      const Point2 b{x(j), y(j)};
      if (b.x == a.x && b.y == a.y) { duplicates++; continue; }
      cell = clip_to_bisector(cell, a, b);
    }
    area[i] = cell.size() >= 3 ? polygon_area(cell) / duplicates : 0.0;
  }
  return area;
}

} // namespace spacc

#endif // SPACC_CORE_GEOMETRY_H

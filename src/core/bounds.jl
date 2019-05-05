
# Bounds #
"Axis aligned bounding box in interval representation"
struct Bounds{N, T <: Point}
  pmin::T
  pmax::T
end

"Initialize a Bounds3 to enclose a single point:" 
Bounds(p::Point) = Bounds(p, p)

"Smallest bounds encapsulating `p1` and `p2`"
Bounds(p1, p2) = Bounds(min.(p1, p2), max.(p1, p2))

"coordinates of `i`th eight corners of the bounding box `b`"
corner(b::Bounds, i) = # fox,e

union(b::Bounds, p::Point) = Bounds(min.(b.pmin, p), max.(b.pmax, p))

union(b1::Bounds, b2::Bounds) = Bounds(min.(b1.pmin, b2.pmin), max.(b1.pmax, b2.pmax))

intersect(b1::Bounds, b2::Bounds) = error("todo")

"Vector along the box diagonal from the minimum point to the maximum point."
diagonal(b::Bounds) = b.pmax - b.pmin

function surfacearea(b::Bounds)
  d = diagonal(b)
  2 * (d.x * d.y + d.x * d.z + d.y * d.z)
end

volume(b::Bounds) = (d = diagonal(b); d.x * d.y * d.z) 

"Index of which of the three axes is longest. This is useful, for example, when deciding which axis to subdivide when building some of the ray-intersection acceleration structures."

"The index of which of the three axes is longest."
function maximumextent(b::Bounds)
  d = diagonal(b)
  if d.x > d.y && d.x > d.z
    1
  elseif d.y > d.z
    2
  else
    3
  end
end

"linearly interpolates between the corners of the box by the given amount in each dimension."
lerp(b::Bounds, t::Point) = Point(lerp(t.x, b.pmin.x, b.pmax.x),
                                  lerp(t.x, b.pmin.x, b.pmax.x),
                                  lerp(t.x, b.pmin.x, b.pmax.x))

"Continuous position of a point relative to the corners of the box, where a point at the minimum corner has offset , a point at the maximum corner has offset , and so forth."
offset(::Bounds, ::Point) = error("todo")

"Sphere that bounds `b`"
boundingsphere(::Bounds) = error("todo")



# template <typename T>
# Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
#     Bounds3<T> ret;
#     ret.pMin = Min(b.pMin, p);
#     ret.pMax = Max(b.pMax, p);
#     return ret;
# }

# template <typename T>
# Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
#     Bounds3<T> ret;
#     ret.pMin = Min(b1.pMin, b2.pMin);
#     ret.pMax = Max(b1.pMax, b2.pMax);
#     return ret;
# }

# template <typename T>
# Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
#     // Important: assign to pMin/pMax directly and don't run the Bounds2()
#     // constructor, since it takes min/max of the points passed to it.  In
#     // turn, that breaks returning an invalid bound for the case where we
#     // intersect non-overlapping bounds (as we'd like to happen).
#     Bounds3<T> ret;
#     ret.pMin = Max(b1.pMin, b2.pMin);
#     ret.pMax = Min(b1.pMax, b2.pMax);
#     return ret;
# }

# template <typename T>
# bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
#     bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
#     bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
#     bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
#     return (x && y && z);
# }

# template <typename T>
# bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
#     return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
#             p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
# }

# template <typename T>
# bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
#     return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
#             p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z);
# }

# template <typename T, typename U>
# inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
#     return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
#                       b.pMax + Vector3<T>(delta, delta, delta));
# }

# distance from point to box; returns zero if point is
# // inside.
# template <typename T, typename U>
# inline Float DistanceSquared(const Point3<T> &p, const Bounds3<U> &b) {
#     Float dx = std::max({Float(0), b.pMin.x - p.x, p.x - b.pMax.x});
#     Float dy = std::max({Float(0), b.pMin.y - p.y, p.y - b.pMax.y});
#     Float dz = std::max({Float(0), b.pMin.z - p.z, p.z - b.pMax.z});
#     return dx * dx + dy * dy + dz * dz;
# }

# template <typename T, typename U>
# inline Float Distance(const Point3<T> &p, const Bounds3<U> &b) {
#     return std::sqrt(DistanceSquared(p, b));
# }

# inline Bounds2iIterator begin(const Bounds2i &b) {
#     return Bounds2iIterator(b, b.pMin);
# }

# inline Bounds2iIterator end(const Bounds2i &b) {
#     // Normally, the ending point is at the minimum x value and one past
#     // the last valid y value.
#     Point2i pEnd(b.pMin.x, b.pMax.y);
#     // However, if the bounds are degenerate, override the end point to
#     // equal the start point so that any attempt to iterate over the bounds
#     // exits out immediately.
#     if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
#         pEnd = b.pMin;
#     return Bounds2iIterator(b, pEnd);
# }

# template <typename T>
# Bounds2<T> Union(const Bounds2<T> &b, const Point2<T> &p) {
#     Bounds2<T> ret;
#     ret.pMin = Min(b.pMin, p);
#     ret.pMax = Max(b.pMax, p);
#     return ret;
# }

# template <typename T>
# Bounds2<T> Union(const Bounds2<T> &b, const Bounds2<T> &b2) {
#     Bounds2<T> ret;
#     ret.pMin = Min(b.pMin, b2.pMin);
#     ret.pMax = Max(b.pMax, b2.pMax);
#     return ret;
# }

# template <typename T>
# Bounds2<T> Intersect(const Bounds2<T> &b1, const Bounds2<T> &b2) {
#     // Important: assign to pMin/pMax directly and don't run the Bounds2()
#     // constructor, since it takes min/max of the points passed to it.  In
#     // turn, that breaks returning an invalid bound for the case where we
#     // intersect non-overlapping bounds (as we'd like to happen).
#     Bounds2<T> ret;
#     ret.pMin = Max(b1.pMin, b2.pMin);
#     ret.pMax = Min(b1.pMax, b2.pMax);
#     return ret;
# }

# template <typename T>
# bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb) {
#     bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
#     bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
#     return (x && y);
# }

# template <typename T>
# bool Inside(const Point2<T> &pt, const Bounds2<T> &b) {
#     return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
#             pt.y <= b.pMax.y);
# }

# template <typename T>
# bool InsideExclusive(const Point2<T> &pt, const Bounds2<T> &b) {
#     return (pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
#             pt.y < b.pMax.y);
# }

# template <typename T, typename U>
# Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
#     return Bounds2<T>(b.pMin - Vector2<T>(delta, delta),
#                       b.pMax + Vector2<T>(delta, delta));
# }

# template <typename T>
# inline bool Bounds3<T>::IntersectP(const Ray &ray, Float *hitt0,
#                                    Float *hitt1) const {
#     Float t0 = 0, t1 = ray.tMax;
#     for (int i = 0; i < 3; ++i) {
#         // Update interval for _i_th bounding box slab
#         Float invRayDir = 1 / ray.d[i];
#         Float tNear = (pMin[i] - ray.o[i]) * invRayDir;
#         Float tFar = (pMax[i] - ray.o[i]) * invRayDir;

#         // Update parametric interval from slab intersection $t$ values
#         if (tNear > tFar) std::swap(tNear, tFar);

#         // Update _tFar_ to ensure robust ray--bounds intersection
#         tFar *= 1 + 2 * gamma(3);
#         t0 = tNear > t0 ? tNear : t0;
#         t1 = tFar < t1 ? tFar : t1;
#         if (t0 > t1) return false;
#     }
#     if (hitt0) *hitt0 = t0;
#     if (hitt1) *hitt1 = t1;
#     return true;
# }

# template <typename T>
# inline bool Bounds3<T>::IntersectP(const Ray &ray, const Vector3f &invDir,
#                                    const int dirIsNeg[3]) const {
#     const Bounds3f &bounds = *this;
#     // Check for ray intersection against $x$ and $y$ slabs
#     Float tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
#     Float tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
#     Float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
#     Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

#     // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
#     tMax *= 1 + 2 * gamma(3);
#     tyMax *= 1 + 2 * gamma(3);
#     if (tMin > tyMax || tyMin > tMax) return false;
#     if (tyMin > tMin) tMin = tyMin;
#     if (tyMax < tMax) tMax = tyMax;

#     // Check for ray intersection against $z$ slab
#     Float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
#     Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

#     // Update _tzMax_ to ensure robust bounds intersection
#     tzMax *= 1 + 2 * gamma(3);
#     if (tMin > tzMax || tzMin > tMax) return false;
#     if (tzMin > tMin) tMin = tzMin;
#     if (tzMax < tMax) tMax = tzMax;
#     return (tMin < ray.tMax) && (tMax > 0);
# }
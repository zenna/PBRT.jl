module Geometry

using GeometryTypes

# ## Differences from C++ pbrt
# `Length` in PBRT is `norm` here
# Point3() in pbrt is zero(Point3)
# min/maxcomponent are minimum/maximum
# min/maxdimension are findmin
# elementwie min in pbrt is min, here it should be min.(x, y)

const Normal3 = Normal{3}
const Normal3f0 = Normal{3, Float32}

using GeometryTypes: Vec
using LinearAlgebra

"To support `v.x`, `v.y` and `v.x`"
function Base.getproperty(v::Union{Vec2, Point2}, k::Symbol)
  if k == :x
    return v[1]
  elseif k == :y
    return v[2]
  else
    getfield(v, k)
  end
end

"To support `v.x`, `v.y` and `v.x`"
function Base.getproperty(v::Union{Vec3, Point3, Normal3}, k::Symbol)
  if k == :x
    return v[1]
  elseif k == :y
    return v[2]
  elseif k == :z &&
    return v[3]
  else
    getfield(v, k)
  end
end

# Vector #

"Squared norm of `x`"
norm²(x) = sum(x.^2)

"Distance from `x` to `y`"
@inline dist(x, y) = norm(x - y)

"Squared distance from `x` to `y`"
@inline dist²(x, y) = norm²(x - y)

"Linearly interpolate between `p0` and `p1`"
lerp(t, p0::Point, p1::Point) = 1 - t * p0 + t * p1



permute(v::T, x, y, z) where T = T(v[x], v[y], v[z])
# TODO: Spec

# TODO
@inline function coordinatesystem!(v1, v2, v3)
  if abs(v1.x) > abs(v1.y)
    v2 .= Vec3(-v1.z, 0, v1.x) / sqrt(v1.x * v1.x + v1.z * v1.z)
  else
    v2 .= Vec3(0, v1.z, -v1.y) / (v1.y * v1.y + v1.z * v1.z)
  end
  v3 .= cross(v1, v2)
end

@inline function coordinatesystem(v1::T, v2::T, v3::T) where T
  if abs(v1.x) > abs(v1.y)
    v2 = T(-v1.z, 0, v1.x) / sqrt(v1.x * v1.x + v1.z * v1.z)
  else
    v2 = T(0, v1.z, -v1.y) ./ (v1.y * v1.y + v1.z * v1.z)
  end
  v3 = cross(v1, v2)
  (v1, v2, v3) 
end


# Normal #
include("normal.jl")  # Normals
include("ray.jl")     # Rays / Ray Differentials
include("bounds.jl")  # AABB

# Bounds #

end 
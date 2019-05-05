"Spheres and PartialSpheres"
abstract type SphereLike{T} end

struct Sphere{T} <: SphereLike{T}
  radius::T
end

area(s::Sphere{T}) where T = 4 * T(pi) * s.radius^2

function objectbound(s::Sphere{T}) where T
  Bounds{T}(Point3{T}(-s.radius), Point3{T}(s.radius))
end

struct PartialSphere{T} <: SphereLike{T}
  radius::T
  zmin::T
  zmax::T
  θmin::T
  θmax::T
  ϕmax::T
end

function objectbound(s::PartialSphere{T}) where T
  Bounds{T}(Point3{T}(-s.radius, -s.radius, s.zmin),
            Point3{T}(s.radius, s.radius, s.zmax))
end

area(s::PartialSphere) = s.ϕmax * s.radius * (s.zmax - s.zmin)

function intersect(r, thit, isect)
  
end


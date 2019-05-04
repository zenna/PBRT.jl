
# Ray #

"Semi-infinite line specified by its origin `o` and direction `d`"
struct Ray{OT <: Point, DT <: Vec, TMAX}
  o::OT
  d::DT
  tmax::TMAX
end

Base.isnan(r::Ray) = isnan(r.o) || isnan(r.d)

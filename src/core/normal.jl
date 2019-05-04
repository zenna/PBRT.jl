const NV3{T} = Union{Normal3{T}, Vec3{T}}

"flip a surface normal `n` so that it lies in the same hemisphere as a given vector `v`"
faceforward(n::NV3{T}, v::NV3{T}) where T = dot(n, v) < zero(T) ? -n : n

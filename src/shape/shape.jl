module Shape

"Abstract type to represent geometric shapes"
abstract type Shape end


function objectbound end
function worldbound end
function intersect end
function intersectp end
function area end
function interaction end
function pdf end
function solidangle end

include("sphere.jl")

end
"""
    ESh(p1::Float64,p2::Float64)

Produces the effect size `h` for two proportions. `p1` and `p2`
are two proportions. `h` is 2*asin(sqrt(p1))-2*asin(sqrt(p2)).
"""
function ESh(p1::Float64,p2::Float64)
    return 2*asin(sqrt(p1)) - 2*asin(sqrt(p2))
end

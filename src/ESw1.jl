"""
    ESw1(P₀::Vector,P₁::Vector)

Compute the effect size `w` for two sets of `k` (number of cells) probabilities.
`P₀` is the vector of `k` proportions under the null hypothesis
and `P₁` is the vector of `k` proportions under the alternative hypothesis.
"""
function ESw1(P0::Vector,P1::Vector)
    return sqrt(sum((P1-P0).^2./P0))
end

"""
    ESw2(P::Array)

Compute the effect size `w`  for a two-way probability table
corresponding to the alternative hypothesis in the chi-squared test
of association in two-way contingency tables
"""
function ESw2(P::Array)
    pcol = sum(P,1)
    prow = sum(P,2)
    P0 = prow * pcol
    return sqrt(sum((P-P0).^2./P0))
end

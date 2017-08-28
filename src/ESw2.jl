function ESw2(P::Array)
    pcol = sum(P,1)
    prow = sum(P,2)
    P0 = prow * pcol
    return sqrt(sum((P-P0).^2./P0))
end

function ESw1(P0::Vector,P1::Vector)
    return sqrt(sum((P1-P0).^2./P0))
end

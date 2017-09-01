# power analysis functions for F2

# u = degree of freedom
# v = degree of freedom
# f2 = effect size
# alpha = Type I error
# power = 1 - Type II error
function powerF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Float64 = 0.0,
    alpha = 0.05)

    if u == 0
        error("`u` must be specified and be greater than 1")
    end

    if v == 0
        error("`v` must be specified and be greater than 1")
    end

    if f2 == 0.0
        error("`f2` must be specified and be positive")
    end

    λ = f2*(u + v +1)
    return ccdf(NoncentralF(u,v,λ),cquantile(FDist(u,v),alpha))
end

#powerF2Test(u=12,v=23,f2=.3)

function samplesizeF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8)

    if alpha == 0.0 || power == 0.0
        error("`alpha` or `power` cannot be zero")
    end

    if u > 0 && v > 0
        error("Either `u` or `v` must be zero")
    end

    if u == 0
        u = v
        u = 0
    end

    return ceil(Int64,fzero(x->powerF2Test(u = u, v = x, f2 = f2, alpha = alpha) - power, 2.0, 10.0^7))
end

function effectsizeF2Test(;
    u::Real = 0,
    v::Real = 0,
    alpha = 0.05,
    power = 0.8)

    return fzero(x->powerF2Test(u = u, v = v, f2 = x, alpha = alpha) - power, 1e-07, 1e+07)
end

function alphaF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2 = 0.0,
    power = 0.8)

    if u < 1 || v < 1
        error("Degress of freedom `u` and `v` must be greater than 1")
    end

    return fzero(x->powerTTest(u = u, v = v, f2 = f2, alpha = x) - power, 1e-10, 1 - 1e-10)
end

function F2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8)

    if sum([x == 0 for x in (u,v,f2,alpha,power)]) != 1
        error("exactly one of u, v, f2, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerF2Test(u = u, v = v, f2 = f2, alpha = alpha)
    elseif alpha == 0.0
        alpha = alphaF2Test(u = u, v = v, f2 = f2, power = power)
    elseif f2 == 0.0
        f2 = effectsizeF2Test(u = u, v = v, alpha = alpha, power = power)
    elseif u == 0
        u = samplesizeF2Test(u = u, v = v, f2 = f2, alpha = alpha, power = power)
    elseif v == 0
        v = samplesizeF2Test(u = u, v = v, f2 = f2, alpha = alpha, power = power)
    end

    return htest(
        string("Multiple regression power calculation"),
        OrderedDict(
            "u" => u,
            "v" => v,
            "f2" => f2,
            "alpha" => alpha,
            "power" => power)
        )
end

# print(F2Test(u=12,v=99,f2=.3,power=0.0))

# power analysis functions for F2

"""
    powerF2Test(u::Real = 0, v::Real = 0, f2::Float64 = 0.0, alpha::Float64 = 0.05)

Compute power for the general linear model for a sample with `u` degress of
freedom for numerator and `v` degrees of freedom for denominator to test the
effect size `f2` at type I error = `alpha` (default: 0.05).
"""
function powerF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Real = 0.0,
    alpha::Real = 0.05)

    check_args(u=u, v=v, f2=f2, alpha=alpha)

    λ = f2*(u + v +1)
    return ccdf(NoncentralF(u,v,λ),cquantile(FDist(u,v),alpha))
end

"""
    samplesizeF2Test(u::Real = 0, v::Real = 0, f2::Float64 = 0.0, alpha::Real = 0.05)

Compute power for the general linear model for a sample with `u` degress of
freedom for numerator and `v` degrees of freedom for denominator to test the
effect size `f2` at type I error = `alpha` (default: 0.05).
"""
function samplesizeF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0)

    check_args(u=u, v=v, f2=f2, alpha=alpha, power=power)

    if u == 0
        return ceil(Int64,fzero(x->powerF2Test(u = x, v = v, f2 = f2, alpha = alpha) - power, 2.0, 10.0^7))
    end
    return ceil(Int64,fzero(x->powerF2Test(u = u, v = x, f2 = f2, alpha = alpha) - power, 2.0, 10.0^7))
end

"""
    effectsizeF2Test(u::Real = 0, v::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8)

Compute the effect size for the general linear model that a sample with `u` degress of
freedom for numerator and `v` degrees of freedom for denominator can test with
power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
"""
function effectsizeF2Test(;
    u::Real = 0,
    v::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0)

    check_args(u=u,v=v,alpha=alpha,power=power)

    return fzero(x->powerF2Test(u = u, v = v, f2 = x, alpha = alpha) - power, 1e-07, 1e+07)
end

"""
    alphaF2Test(u::Real = 0, v::Real = 0, f2::Float64 = 0.0, power::Float64 = 0.8)

Compute the probability of type I error that the general linear model with
a sample with `u` degress of freedom for numerator and `v` degrees of freedom
for denominator has in testin the effect size `f2` with power > `power` (default: 0.8).
"""
function alphaF2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Real = 0.0,
    power::Real = 0.0)

    check_args(u=u, v=v, f2=f2, power=power)

    return fzero(x->powerTTest(u = u, v = v, f2 = f2, alpha = x) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.F2Test(u::Real = 0, v::Real = 0, f2::Float64 = 0.0, alpha::Float64 = 0.0, power::Float64 = 0.0)

Compute one of the test parameters such as sample size with `u` degrees of freedom
for the numerator and `v` degrees of freedom for the denominator,
effect size (`f2`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
"""
function F2Test(;
    u::Real = 0,
    v::Real = 0,
    f2::Real = 0.0,
    alpha::Real = 0.0,
    power::Real = 0.0)

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

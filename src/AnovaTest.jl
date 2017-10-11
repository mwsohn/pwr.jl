# power analysis functions for ANOVA

# k = number of groups
# n = number of observations per group
# f = effect size
# alpha = Type I error
# power = 1 - Type II error

"""
    powerAnovaTest(k::Real = 0, n::Real = 0, f::Real = 0.0, alpha::Real = 0.05)

Estimates power of a sample with `n` observations in each of the `k` groups
to detect an effect size `f` at type I error = `alpha` (default: 0.05).
"""
function powerAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f::Real = 0.0,
    alpha::Real = 0.05)

    check_args(k=k,n=n,f=f,alpha=alpha)

    λ = k*n*f^2
    return ccdf(NoncentralF(k-1,(n-1)*k,λ),cquantile(FDist(k-1,(n-1)*k),alpha))
end

"""
    samplesizeAnovaTest(k::Real = 0, f::Real = 0.0, alpha::Real = 0.05, power::Real = 0.8)

Estimates the minimum sample size required to achieve power > `power` (default: 0.8) in
each of the `k` groups to detect an effect size `f` at type I error = `alpha` (default: 0.05).
"""
function samplesizeAnovaTest(;
    k::Real = 0,
    f::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.8)

    check_args(k=k,f=f,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerAnovaTest(k = k, n = x, f = f, alpha = alpha) - power, 2.0, 10.0^7))
end

"""
    groupsizeAnovaTest(n::Real = 0, f::Real = 0.0, alpha::Real = 0.05, power::Real = 0.8)

Estimates the number of groups required to achieve power > `power` (default: 0.8) with
`n` observations in each group to detect an effect size `f` at type I error = `alpha` (default: 0.05).
"""
function groupsizeAnovaTest(;
    n::Real = 0,
    f::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.8)

    check_args(n=n,f=f,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerAnovaTest(k = k, n = x, f = f, alpha = alpha) - power, 2.0, 10.0^7))
end

"""
    effectsizeAnovaTest(k::Real = 0, n::Real = 0.0, alpha::Real = 0.05, power::Real = 0.8)

Estimates the effect size that a sample with `n` observations each in `k` groups
can test at power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
"""
function effectsizeAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.8)

    check_args(k=k,n=n,alpha=alpha,power=power)

    return fzero(x->powerAnovaTest(k = k, n = n, f = x, alpha = alpha) - power, 1e-07, 1e+07)
end


"""
    alphaAnovaTest(k::Real = 0, n::Real = 0.0, f::Real = 0.0, power::Real = 0.8)

Estimates the minimum type I error (α) that a sample with `n` observations
each in `k` groups can achieve at power > `power` (default: 0.8) to detect an effect size `f`.
"""
function alphaAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f::Real = 0.0,
    power::Real = 0.8)

    check_args(k=k,f=f,n=n,power=power)

    return fzero(x->powerTTest(n = n, d = d, f = f, alpha = x) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.AnovaTest(k::Real = 0, n::Real = 0.0, f::Real = 0.0, alpha::Real = 0.0, power::Real = 0.0)

Estimates one of the test parameters such as the number of groups (`k`), observations in each group (`n`),
effect size (`f`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
"""
function AnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f::Real = 0.0,
    alpha::Real = 0.0,
    power::Real = 0.0)

    if sum([x == 0 for x in (k,n,f,alpha,power)]) != 1
        error("exactly one of k, n, f, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerAnovaTest(k = k, n = n, f = f, alpha = alpha)
    elseif alpha == 0.0
        alpha = alphaAnovaTest(k = k, n = n, f = f, power = power)
    elseif f == 0.0
        f = effectsizeAnovaTest(k = k, n = n, alpha = alpha, power = power)
    elseif n == 0
        n = samplesizeAnovaTest(k = k, f = f, alpha = alpha, power = power)
    elseif k == 0
        k = groupsizeAnovaTest(n = n, f = f, alpha = alpha, power = power)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")

    note = "`n` is the number in each group"

    return htest(
        string("Balanced one-way analysis of variance power calculation"),
        OrderedDict(
            "k" => k,
            "n" => n,
            "f" => f,
            "alpha" => alpha,
            "power" => power,
            "note" => note)
        )
end

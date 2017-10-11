# power analysis functions for mean for normal distribution with known variance

"""
    powerRTest(r::Real = 0, n::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of a sample of `n` observations to test the effect size `r`
at type I error = `alpha` (default: 0.05).
The effect size `r` is a linear correlation coefficient. The option
`alternative` is `two.sided` for H₀: r = 0 vs H₁: r ≠ 0, `less` for H₁: r < 0, and `greater`
for H₁: r > 0.
"""
function powerRTest(;
    r::Real = 0.0,
    n::Real = 0,
    alpha::Real = 0.05,
    alternative::String = "two")

    check_args(r=r,n=n,alpha=alpha)

    if alternative == "less"
        tside = 1
        r = -r
    elseif alternative in ("two","two.sided","two-sided","two sided")
        tside = 2
        r = abs(r)
    elseif alternative == "greater"
        tside = 1
    end

    if tside == 1
        ttt = cquantile(TDist(n-2),alpha)
        rc = sqrt(ttt^2 / (ttt^2 + n - 2))
        zr = atanh(r) + r / (2 * (n - 1))
        zrc = atanh(rc)
        return cdf(Normal(),(zr - zrc) * sqrt(n-3))
    end
    ttt = cquantile(TDist(n-2),alpha/2)
    rc = sqrt(ttt^2 / (ttt^2 + n - 2))
    zr = atanh(r) + r / (2 * (n - 1))
    zrc = atanh(rc)
    return cdf(Normal(),(zr - zrc) * sqrt(n-3)) + cdf(Normal(),(-zr - zrc) * sqrt(n-3))
end

"""
    samplesizeRTest(r::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size to detect the effect size `r` with power > `power`
at type I error = `alpha` (default: 0.05).
The effect size `r` is a linear correlation coefficient. The option
`alternative` is `two.sided` for H₀: r = 0 vs H₁: r ≠ 0, `less` for H₁: r < 0, and `greater`
for H₁: r > 0.
"""
function samplesizeRTest(;
    r::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two")

    check_args(r=r,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerRTest(n = x, r = r, alpha = alpha, alternative = alternative) - power, 4.0+1e-10, 1e+09))
end

"""
    effectsizeRTest(n::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that a sample size `n` can detect with power > `power`
at type I error = `alpha` (default: 0.05).
The effect size `r` is a linear correlation coefficient. The option
`alternative` is `two.sided` for H₀: r = 0 vs H₁: r ≠ 0, `less` for H₁: r < 0, and `greater`
for H₁: r > 0.
"""
function effectsizeRTest(;
    n::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerRTest(n = n, r = x, alpha = alpha, alternative = alternative) - power,.001,100)
end

"""
    alphaRTest(r::Real = 0, n::Real = 0, power::Float64 = 0.8, alternative = "two.sided")

Compute the probability of type I error that the sample size `n` has in detecting the
effect size `r` with power > `power` at type I error = `alpha` (default: 0.05).
The effect size `r` is a linear correlation coefficient. The option
`alternative` is `two.sided` for H₀: r = 0 vs H₁: r ≠ 0, `less` for H₁: r < 0, and `greater`
for H₁: r > 0.
"""
function alphaRTest(;
    n::Real = 0,
    r::Real = 0.0,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(r=r,n=n,power=power)

    return fzero(x->powerRTest(n = n, r = r, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.RTest(r::Real = 0, n::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n`),
effect size (`r`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size `r` is a linear correlation coefficient. The option
`alternative` is `two.sided` for H₀: r = 0 vs H₁: r ≠ 0, `less` for H₁: r < 0, and `greater`
for H₁: r > 0.
"""
function RTest(;
    n::Real = 0,
    r::Real = 0.0,
    alpha::Float64 = 0.0,
    power::Float64 = 0.0,
    alternative::String = "two")

    if sum([x == 0 for x in (n,r,alpha,power)]) != 1
        error("exactly one of n, r, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerRTest(n = n, r = r, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaRTest(n = n, r = r, power = power, alternative = alternative)
    elseif r == 0.0
        r = effectsizeRTest(n = n, alpha = alpha, power = power, alternative = alternative)
    elseif n == 0
        n = samplesizeRTest(r = r, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided","two.sided" => "two-sided", "two sided" => "two-sided", "less" => "less", "greater" => "greater")

    return htest(
        string("Approximate correlation power calculation (arctangh transformation)"),
        OrderedDict(
            "n" => n,
            "r" => r,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative])
        )
end

# power analysis functions for mean for normal distribution with known variance

"""
    powerNormTest(d::Real = 0, n::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of a sample of `n` observations to test the effect size `d`
at type I error = `alpha` (default: 0.05). The effect size d = μ - μ₀. The option
`alternative` is `two.sided` for H₀: μ = μ₀ vs H₁: μ ≠ μ₀, `less` for H₁: μ < μ₀, and `greater`
for H₁: μ > μ₀.
"""
function powerNormTest(;
    d::Real = 0.0,
    n::Real = 0,
    alpha::Float64 = 0.05,
    alternative::String = "two")

    check_args(d=d,n=n,alpha=alpha)

    if alternative == "less"
        tside = 1
    elseif alternative in ("two","two.sided","two-sided")
        tside = 2
        d = abs(d)
    elseif alternative == "greater"
        tside = 3
    end

    if tside == 1
        return cdf(Normal(),quantile(Normal(),alpha) - d*sqrt(n))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - d*sqrt(n)) + cdf(Normal(),quantile(Normal(),alpha/2) - d*sqrt(n))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha) - d*sqrt(n))
    end
end

"""
    samplesizeNormTest(d::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size that can test the effect size `d` with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05). The effect size d = μ - μ₀. The option
`alternative` is `two.sided` for H₀: μ = μ₀ vs H₁: μ ≠ μ₀, `less` for H₁: μ < μ₀, and `greater`
for H₁: μ > μ₀.
"""
function samplesizeNormTest(;
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two")

    check_args(d=d,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerNormTest(n = x, d = d, alpha = alpha, alternative = alternative) - power, 2.0, 10.0^7))
end

"""
    effectsizeNormTest(alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that a sample size `n` can test with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05). The option
`alternative` is `two.sided` for H₀: μ = μ₀ vs H₁: μ ≠ μ₀, `less` for H₁: μ < μ₀, and `greater`
for H₁: μ > μ₀.
"""
function effectsizeNormTest(;
    n::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerNormTest(n = n, d = x, alpha = alpha, alternative = alternative) - power,.001,100)
end

"""
    alphaNormTest(d::Real = 0, n::Real = 0, power::Float64 = 0.8, alternative = "two.sided")

Compute the probability of type I error that a sample size `n` has in testing
the effect size `d` with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
The effect size d = μ - μ₀. The option `alternative` is `two.sided` for H₀: μ = μ₀ vs H₁: μ ≠ μ₀,
`less` for H₁: μ < μ₀, and `greater` for H₁: μ > μ₀.
"""
function alphaNormTest(;
    n::Real = 0,
    d::Real = 0.0,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(d=d,n=n,power=power)

    return fzero(x->powerNormTest(n = n, d = d, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.NormTest(d::Real = 0, n::Real = 0, alpha::Float64 = 0.0, power::Float64 = 0.0, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n`),
effect size (`d`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size d = μ - μ₀. The option `alternative` is `two.sided` for H₀: μ = μ₀ vs H₁: μ ≠ μ₀,
`less` for H₁: μ < μ₀, and `greater` for H₁: μ > μ₀.
"""
function NormTest(;
    n::Real = 0,
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two")

    if sum([x == 0 for x in (n,d,alpha,power)]) != 1
        error("exactly one of n, d, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerNormTest(n = n, d = d, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaNormTest(n = n, d = d, power = power, alternative = alternative)
    elseif d == 0.0
        d = effectsizeNormTest(n = n, alpha = alpha, power = power, alternative = alternative)
    elseif n == 0
        n = samplesizeNormTest(d = d, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")

    return htest(
        string("Mean power calculation for normal distribution with known variance"),
        OrderedDict(
            "n" => n,
            "d" => d,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative]
            )
        )
end

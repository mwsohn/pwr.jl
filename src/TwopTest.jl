# power analysis functions for difference in two proportions

"""
    powerTwoPTest(h::Real = 0, n::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of two equal size samples of `n` observations to test the effect size `h`
for a difference between two independent proportions (P) at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function powerTwoPTest(;
    h::Real = 0,
    n::Real = 0,
    alpha::Float64 = 0.05,
    alternative::String = "two"
    )

    check_args(h=h,n=n,alpha=alpha)

    if alternative == "less"
        tside = 1
    elseif alternative in ("two","two.sided","two-sided")
        tside = 2
        h = abs(h)
    elseif alternative == "greater"
        tside = 3
    else
        error("`alternative` must be `less`, `two`, or `greater`.")
    end

    if tside == 1
        return cdf(Normal(),quantile(Normal(),alpha) - h*sqrt(n/2))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - h*sqrt(n/2)) + cdf(Normal(),quantile(Normal(),alpha/2) - h*sqrt(n/2))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha) - h*sqrt(n/2))
    end
    error("internal error")
end

"""
    samplesizeTwoPTest(h::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size `n` of two equal size samples to test the effect size `h`
for a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function samplesizeTwoPTest(;
    h::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(h=h,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTwoPTest(h = h, n = x, alpha = alpha, alternative = alternative) - power, 2 + 1e-10, 1e+09))
end

"""
    effectsizeTwoPTest(n::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that two equal size samples of `n` observations can detect for
a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function effectsizeTwoPTest(;
    n::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    if alternative == "less"
        return fzero(x->powerTwoPTest(h = x, n = n, alpha = alpha, alternative = alternative) - power, -10.0, 5.0)
    elseif alternative == "two"
        return fzero(x->powerTwoPTest(h = x, n = n, alpha = alpha, alternative = alternative) - power, 1e-10, 10.0)
    elseif side == "greater"
        return fzero(x->powerTwoPTest(h = x, n = n, alpha = alpha, alternative = alternative) - power, -5.0, 10.0)
    end
end

"""
    alphaTwoPTest(h::Real = 0, n::Real = 0, power::Float64 = 0.8, alternative = "two.sided")

Compute the probability of type I error that two equal size samples of `n` observations have in detecting
the effect size `h` for a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function alphaTwoPTest(;
    h::Real = 0,
    n::Real = 0,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(h=h,n=n,power=power)

    return fzero(x->powerTwoPTest(h = h, n = n, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.TwoPTest(h::Real = 0, n::Real = 0, alpha::Float64 = 0.0, power::Float64 = 0.0, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n`),
effect size (`h`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)).
The option `alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function TwoPTest(;
    h::Real = 0,
    n::Real = 0,
    alpha::Float64 = 0.0,
    power::Float64 = 0.0,
    alternative::String = "two"
    )
    if sum([x == 0 for x in (h,n,alpha,power)]) != 1
        error("exactly one of `h`, `n`, `power`, and `alpha` must be zero")
    end

    if power == 0.0
        power = powerTwoPTest(h = h, n = n, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaTwoPTest(h = h, n = n, power = power, alternative = alternative)
    elseif h == 0
        h = effectsizeTwoPTest(n = n, alpha = alpha, power = power, alternative = alternative)
    elseif n == 0
        n = samplesizeTwoPTest(h = h, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided","less" => "less", "greater" => "greater")

    return htest(
        "Difference of proportion power calculation for binomial distribution (same sample)",
        OrderedDict(
            "h" => h,
            "n" => n,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative])
        )
end

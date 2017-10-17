# power analysis functions for difference in two proportions in two different sample sizes

"""
    powerTwoP2nTest(h::Real = 0, n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of two samples of `n1` and `n2` observations (n1 ≠ n2) to test the effect size `h`
for a difference between two independent proportions (P) at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function powerTwoP2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    alternative::String = "two"
    )

    check_args(h=h,n1=n1,n2=n2,alpha=alpha)

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
        return cdf(Normal(),quantile(Normal(),alpha) - h*sqrt(n1*n2/(n1 + n2)))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - h*sqrt(n1*n2/(n1 + n2))) + cdf(Normal(),quantile(Normal(),alpha/2) - h*sqrt(n1*n2/(n1 + n2)))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha) - h*sqrt(n1*n2/(n1 + n2)))
    end
    error("internal error")
end

"""
    samplesizeTwoP2nTest(h::Real = 0, n1::Real = 0,alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size `n2` of two unequal size samples to test the effect size `h`
for a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function samplesizeTwoP2nTest(;
    h::Real = 0,
    n1::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(h=h,n1=n1,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTwoP2nTest(h = h, n1 = n1, n2 = x, alpha = alpha, alternative = alternative) - power, 2 + 1e-10, 1e+09))
end

"""
    effectsizeTwoP2nTest(n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that two equal size samples of `n` observations can detect for
a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function effectsizeTwoP2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(n1=n1,n2=n2,alpha=alpha,power=power)

    return fzero(x->powerTwoP2nTest(h = x, n1 = n1, n2 = n2, alpha = alpha, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    alphaTwoP2nTest(h::Real = 0, n1::Real = 0, n2::Real = 0, power::Float64 = 0.05, alternative = "two.sided")

Compute the probability of type I error that two unequal size samples of `n1` and `n2` observations have in detecting
the effect size `h` for a difference between two independent proportions (P) with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)). The option
`alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function alphaTwoP2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    power::Float64 = 0.0,
    alternative::String = "two"
    )

    check_args(h=h,n1=n1,n2=n2,power=power)

    return fzero(x->powerTwoP2nTest(h = h, n1 = n1, n2 = n2, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    TwoP2nTest(h::Real = 0, n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.0, power::Float64 = 0.0, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n2`),
effect size (`h`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size `h` = ϕ₁ - ϕ₂ where ϕ = 2*asin(sqrt(P)).
The option `alternative` is `two.sided` for H₀: ϕ₁ = ϕ₂ vs H₁: ϕ₁ ≠ ϕ₂, `less` for H₁: ϕ₁ < ϕ₂, and `greater`
for H₁: ϕ₁ > ϕ₂.
"""
function TwoP2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.0,
    alternative::String = "two"
    )
    if sum([x == 0 for x in (h,n1,n2,alpha,power)]) != 1
        error("exactly one of `h`, `n2`, `power`, and `alpha` must be zero")
    end

    if n1 == 0
        error("`n1` cannot be zero")
    end

    if power == 0.0
        power = powerTwoP2nTest(h = h, n1 = n1, n2 = n2, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaTwoP2nTest(h = h, n1 = n1, n2 = n2, power = power, alternative = alternative)
    elseif h == 0
        h = effectsizeTwoP2nTest(n1 = n1, n2 = n2, alpha = alpha, power = power, alternative = alternative)
    elseif n2 == 0
        n2 = samplesizeTwoP2nTest(h = h, n1 = n1, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")

    note = "different sample sizes"

    return htest(
        "Difference of proportion power calculation for binomial distribution (different sample)",
        OrderedDict(
            "h" => h,
            "n1" => n1,
            "n2" => n2,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative],
            "note" => note)
        )
end

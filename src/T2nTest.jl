# power analysis functions for two sample t-test

"""
    powerT2nTest(d::Real = 0, n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of an unbalanced two samples of `n1` and `n2` observations to test the effect size `d`
at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂, `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest, `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function powerT2nTest(;
    d::Real = 0.0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    alternative::String = "two")

    check_args(d=d,n1=n1,n2=n2,alpha=alpha)

    if alternative == "less"
        tside = ttside = 1
    elseif alternative in ("two","twosided","two.sided","two-sided","two sided")
        tside = ttside = 2
        d = abs(d)
    elseif alternative == "greater"
        ttside = 3
        tside = 1
    end

    ν = ceil(Int64,n1 + n2 - 2)
    λ = d*(1/sqrt(1/n1 + 1/n2))
    if ttside == 1
        return cdf(NoncentralT(ν,λ),quantile(TDist(ν),alpha))
    elseif ttside == 2
        qu = cquantile(TDist(ν),alpha/tside)
        return ccdf(NoncentralT(ν,λ),qu) + cdf(NoncentralT(ν,λ),-qu)
    elseif ttside == 3
        return ccdf(NoncentralT(ν,λ),cquantile(TDist(ν),alpha))
    end
end

"""
    samplesizeT2nTest(d::Real = 0, n1::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size `n2` for a unbalanced two samples of which one has `n1` observations
to test the effect size `d` with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function samplesizeT2nTest(;
    n1::Real = 0,
    d::Real = 0.0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two")

    check_args(d=d,n1=n1,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerT2nTest(n1 = n1, n2 = x, d = d, alpha = alpha, alternative = alternative) - power, 2+1e-09, 1e+10))
end

"""
    effectsizeT2nTest(n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that unbalanced two samples of `n1` and `n2` observations can detect
with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function effectsizeT2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(n1=n1,alpha=alpha,power=power)

    return fzero(x -> powerT2nTest(n1 = n1, n2 = n2, d = x, alpha = alpha, alternative = alternative) - power,.001,100)
end


"""
    alphaT2nTest(n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the probability of type I error that unbalanced two samples of `n1` and `n2` observations have
in detecting the effect size `d` with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function alphaT2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    d::Real = 0.0,
    power::Float64 = 0.0,
    alternative::String = "two"
    )

    check_args(d=d,n1=n1,n2=n2,power=power)

    return fzero(x->powerT2nTest(n1 = n1, n2 = n2, d = d, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end


"""
    pwr.T2nTest(d::Real = 0, n1::Real = 0, n2::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n2`),
effect size (`d`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function T2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    d::Real = 0.0,
    alpha::Float64 = 0.0,
    power::Float64 = 0.0,
    alternative::String = "two")

    if sum([x == 0 for x in (n1,n2,d,alpha,power)]) != 1
        error("exactly one of `n2`, `d`, `power`, and `alpha` must be zero")
    end

    if power == 0.0
        power = powerT2nTest(n1 = n1, n2 = n2, d = d, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaT2nTest(n1 = n1, n2 = n2, d = d, power = power, alternative = alternative)
    elseif d == 0.0
        d = effectsizeT2nTest(n1 = n1, n2 = n2, alpha = alpha, power = power, alternative = alternative)
    elseif n2 == 0
        n2 = samplesizeT2nTest(n1 = n1, d = d, alpha = alpha, power = power, alternative = alternative)
    end

    stype = Dict("onesample" => "One-sample", "one.sample" => "One-sample", "one sample" => "One-sample",
        "twosample" => "Two-sample", "two.sample" => "Two-sample", "two sample" => "Two-sample",
        "paired" => "Paired")
    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")
    note = ""

    return htest(
        string("T-test power calculation"),
        OrderedDict(
            "n1" => n1,
            "n2" => n2,
            "d" => d,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative]
            )
        )
end

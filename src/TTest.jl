# power analysis functions for t-test

"""
    powerTTest(d::Real = 0, n::Real = 0, alpha::Float64 = 0.05, alternative = "two.sided")

Compute power of a sample of `n` observations to test the effect size `d`
at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂, `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest, `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function powerTTest(;
    d::Real = 0.0,
    n::Real = 0,
    alpha::Float64 = 0.05,
    sampletype::String = "onesample",
    alternative::String = "two")

    check_args(d=d,n=n,alpha=alpha)

    if lowercase(sampletype) in ("onesample","one.sample","one-sample","one sample","paired")
        tsample = 1
    elseif lowercase(sampletype) in ("twosample","two.sample","two-sample","two sample")
        tsample = 2
    end

    if lowercase(alternative) == "less"
        tside = ttside = 1
    elseif lowercase(alternative) in ("two","two.sided","two-sided","two sided")
        tside = ttside = 2
        d = abs(d)
    elseif lowercase(alternative) == "greater"
        ttside = 3
        tside = 1
    end

    ν = ceil(Int64,(n-1)*tsample)
    λ = sqrt(n/tsample)*d
    if ttside == 1
        return cdf(NoncentralT(ν,λ),quantile(TDist(ν),alpha/tside))
    elseif ttside == 2
        qu = cquantile(TDist(ν),alpha/tside)
        return ccdf(NoncentralT(ν,λ),qu) + cdf(NoncentralT(ν,λ),-qu)
    elseif ttside == 3
        return ccdf(NoncentralT(ν,λ),cquantile(TDist(ν),alpha/tside))
    end
end

"""
    samplesizeTTest(d::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the sample size to test the effect size `d` with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function samplesizeTTest(;
    d::Real = 0.0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    sampletype::String = "onesample",
    alternative::String = "two")

    check_args(d=d,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTTest(n = x, d = d, alpha = alpha, sampletype = sampletype, alternative = alternative) - power, 2.0, 10.0^7))
end

"""
    effectsizeTTest(n::Real = 0, alpha::Float64 = 0.05, power::Float64 = 0.8, alternative = "two.sided")

Compute the effect size that a sample of `n` can detect with power > `power` (default: 0.8)
at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function effectsizeTTest(;
    n::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    sampletype::String = "onesample",
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerTTest(n = n, d = x, alpha = alpha, sampletype = sampletype, alternative = alternative) - power,.001,100)
end

"""
    alphaTTest(d::Real = 0, n::Real = 0, power::Float64 = 0.8, alternative = "two.sided")

Compute the probability of type I error that a sample of `n` has in detecting the
effect size `d` with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function alphaTTest(;
    n::Real = 0,
    d::Real = 0.0,
    power::Float64 = 0.8,
    sampletype::String = "onesample",
    alternative::String = "two"
    )

    check_args(d=d,n=n,power=power)

    return fzero(x->powerTTest(n = n, d = d, alpha = x, sampletype = sampletype, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

"""
    pwr.TTest(d::Real = 0, n::Real = 0, alpha::Float64 = 0.0, power::Float64 = 0.0, alternative = "two.sided")

Compute one of the test parameters such as sample size (`n`),
effect size (`d`), type I error (`alpha`), and power (`power`).
The parameter to be estimated must be set to zero (default).
The effect size `d` is the difference in means divided by the pooled standard deviation. The option
`alternative` is `two.sided` for H₀: μ₁ = μ₂ vs H₁: μ₁ ≠ μ₂ (default), `less` for H₁: μ₁ < μ₂, and `greater`
for H₁: μ₁ > μ₂. `sampletype` is `onesample` for one-sample ttest (default), `twosample` for two-sample ttest,
and `paired` for two-sample paired ttest.
"""
function TTest(;
    n::Real = 0,
    d::Real = 0.0,
    alpha::Float64 = 0.0,
    power::Float64 = 0.0,
    sampletype::String = "onesample",
    alternative::String = "two")

    if sum([x == 0 for x in (n,d,alpha,power)]) != 1
        error("exactly one of n, d, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerTTest(n = n, d = d, alpha = alpha, sampletype = sampletype, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaTTest(n = n, d = d, power = power, sampletype = sampletype, alternative = alternative)
    elseif d == 0.0
        d = effectsizeTTest(n = n, alpha = alpha, power = power, sampletype = sampletype, alternative = alternative)
    elseif n == 0
        n = samplesizeTTest(d = d, alpha = alpha, power = power, sampletype = sampletype, alternative = alternative)
    end

    stype = Dict("onesample" => "One-sample","one.sample" => "One-sample","one sample" => "One-sample",
        "one-sample" => "One-sample","twosample" => "Two-sample","two.sample" => "Two-sample",
        "two sample" => "Two-sample","two-sample" => "Two-sample","paired" => "Paired")
    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")

    note = ""
    if alternative == "paired"
        note = "`n` is number of pairs"
    elseif alternative in ("two","two.sided","two-sided","two sided")
        note = "`n` is number in each group"
    end

    return htest(
        string(stype[sampletype]," t-test power calculation"),
        OrderedDict(
            "n" => n,
            "d" => d,
            "alpha" => alpha,
            "power" => power,
            "sampletype" => stype[sampletype],
            "alternative" => alt[alternative],
            "note" => note)
        )
end

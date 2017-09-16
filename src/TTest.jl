# power analysis functions for t-test

function powerTTest(;
    d::Real = 0.0,
    n::Real = 0,
    alpha::Real = 0.05,
    sampletype::String = "onesample",
    alternative::String = "two")

    check_args(d=d,n=n,alpha=alpha)

    if lowercase(sampletype) in ("onesample","one.sample","one-sample","one sample","paired")
        tsample = 1
    elseif lowercase(sampletype) in ("twosample","two.sample","two-sample","two sample")
        tsample = 2
    end

    if alternative == "less"
        tside = ttside = 1
    elseif alternative in ("two","two.sided","two-sided","two sided")
        tside = ttside = 2
        d = abs(d)
    elseif alternative == "greater"
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

function samplesizeTTest(;
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    sampletype::String = "onesample",
    alternative::String = "two")

    check_args(d=d,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTTest(n = x, d = d, alpha = alpha, sampletype = sampletype, alternative = alternative) - power, 2.0, 10.0^7))
end

function effectsizeTTest(;
    n::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    sampletype::String = "onesample",
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerTTest(n = n, d = x, alpha = alpha, sampletype = sampletype, alternative = alternative) - power,.001,100)
end

function alphaTTest(;
    n::Real = 0,
    d::Real = 0.0,
    power::Real = 0.0,
    sampletype::String = "onesample",
    alternative::String = "two"
    )

    check_args(d=d,n=n,power=power)

    return fzero(x->powerTTest(n = n, d = d, alpha = x, sampletype = sampletype, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

function TTest(;
    n::Real = 0,
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
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

    stype = Dict("onesample" => "one-sample","one.sample" => "one-sample","one sample" => "one-sample",
        "twosample" => "two-sample","two.sample" => "two-sample","two sample" => "two-sample",
        "paired" => "paired")
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

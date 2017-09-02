# power analysis functions for t-test

using Distributions, Roots

function powerTTest(;
    d::Float64 = 0.0,
    n::Real = 0,
    alpha = 0.05,
    sampletype::String = "onesample",
    sided::String = "two")

    check_args(d=d,n=n,alpha=alpha)

    if sampletype in ("onesample","paired")
        tsample = 1
    elseif sampletype == "twosample"
        tsample = 2
    end

    if sided == "less"
        tside = ttside = 1
    elseif sided in ("two","two-sided")
        tside = ttside = 2
        d = abs(d)
    elseif sided == "greater"
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
    d::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sampletype::String = "onesample",
    sided::String = "two")

    check_args(d=d,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTTest(n = x, d = d, alpha = alpha, sampletype = sampletype, sided = sided) - power, 2.0, 10.0^7))
end

function effectsizeTTest(;
    n::Int64 = 0,
    alpha = 0.05,
    power = 0.8,
    sampletype::String = "onesample",
    sided::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerTTest(n = n, d = x, alpha = alpha, sampletype = sampletype, sided = sided) - power,.001,100)
end

function alphaTTest(;
    n = 0,
    d = 0.0,
    power = 0.8,
    sampletype::String = "onesample",
    sided::String = "two"
    )

    check_args(d=d,n=n,power=power)

    return fzero(x->powerTTest(n = n, d = d, alpha = x, sampletype = sampletype, sided = sided) - power, 1e-10, 1 - 1e-10)
end

function TTest(;
    n::Real = 0,
    d::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sampletype::String = "onesample",
    sided::String = "two")

    if sum([x == 0 for x in (n,d,alpha,power)]) != 1
        error("exactly one of n, d, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerTTest(n = n, d = d, alpha = alpha, sampletype = sampletype, sided = sided)
    elseif alpha == 0.0
        alpha = alphaTTest(n = n, d = d, power = power, sampletype = sampletype, sided = sided)
    elseif d == 0.0
        d = effectsizeTTest(n = n, alpha = alpha, power = power, sampletype = sampletype, sided = sided)
    elseif n == 0
        n = samplesizeTTest(d = d, alpha = alpha, power = power, sampletype = sampletype, sided = sided)
    end

    stype = Dict("onesample" => "One-sample", "twosample" => "Two-sample", "paired" => "Paired")
    alt = Dict("two" => "two-sided", "less" => "less", "greater" => "greater")

    note = ""
    if sided == "paired"
        note = "`n` is number of pairs"
    elseif sided == "two"
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
            "alternative" => alt[sided],
            "note" => note)
        )
end

#print(TTest(n=0,d=.2))

# power analysis functions for two sample t-test

function powerT2Test(;
    d::Real = 0.0,
    n1::Real = 0,
    n2::Real = 0,
    alpha = 0.05,
    alternative::String = "two")

    check_args(d = d, n1 = n1, n2 = n2, alpha = alpha)

    if alternative == "less"
        tside = ttside = 1
    elseif alternative in ("two","two.sided","two-sided")
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

function samplesizeT2Test(;
    n1::Real = 0,
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two")

    check_args(d = d, n1 = n1, alpha = alpha, power = power)

    return fzero(x->powerT2Test(n1 = n1, n2 = x, d = d, alpha = alpha, alternative = alternative) - power, 2+1e-09, 1e+10)
end

function effectsizeT2Test(;
    n1::Real = 0,
    n2::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(n1 = n1, n2 = n2, alpha = alpha, power = power)

    return fzero(x -> powerT2Test(n1 = n1, n2 = n2, d = x, alpha = alpha, alternative = alternative) - power,.001,100)
end

#effectsizeT2Test(n1=1500,n2=50)

function alphaT2Test(;
    n1::Real = 0,
    n2::Real = 0,
    d::Real = 0.0,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(d = d, n1 = n1, n2 = n2, power = power)

    return fzero(x->powerT2Test(n1 = n1, n2 = n2, d = d, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

function T2Test(;
    n1::Real = 0,
    n2::Real = 0,
    d::Real = 0.0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two")

    if sum([x == 0 for x in (n1,n2,d,alpha,power)]) != 1
        error("exactly one of `n2`, `d`, `power`, and `alpha` must be zero")
    end

    if power == 0.0
        power = powerT2Test(n1 = n1, n2 = n2, d = d, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaT2Test(n1 = n1, n2 = n2, d = d, power = power, alternative = alternative)
    elseif d == 0.0
        d = effectsizeT2Test(n1 = n1, n2 = n2, alpha = alpha, power = power, alternative = alternative)
    elseif n2 == 0
        n2 = samplesizeT2Test(n1 = n1, d = d, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided", "less" => "less", "greater" => "greater")
    note = ""

    return htest(
        "T-test power calculation",
        OrderedDict(
            "n1" => n1,
            "n2" => n2,
            "d" => d,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative],
            "note" => note)
            )
end

#T2Test(n1=200,n2=0,d=.2)
#samplesizeT2Test(n1=100,d=.3)

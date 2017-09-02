# power analysis functions for mean for normal distribution with known variance

function powerNormTest(;
    d::Float64 = 0.0,
    n::Real = 0,
    alpha = 0.05,
    sided::String = "two")

    check_args(d=d,n=n,alpha=alpha)

    if sided == "less"
        tside = 1
    elseif sided in ("two","two-sided")
        tside = 2
        d = abs(d)
    elseif sided == "greater"
        tside = 3
    end

    if tside == 1
        return cdf(Normal(),quantile(Normal(),alpha) - d*sqrt(n))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - d*sqrt(n))
            + cdf(Normal(),quantile(Normal(),alpha/2) - d*sqrt(n))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha) - d*sqrt(n))
    end
end

function samplesizeNormTest(;
    d::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two")

    check_args(d=d,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerNormTest(n = x, d = d, alpha = alpha, sided = sided) - power, 2.0, 10.0^7))
end

function effectsizeNormTest(;
    n::Int64 = 0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerNormTest(n = n, d = x, alpha = alpha, sided = sided) - power,.001,100)
end

function alphaNormTest(;
    n = 0,
    d = 0.0,
    power = 0.8,
    sided::String = "two"
    )

    check_args(d=d,n=n,power=power)

    return fzero(x->powerNormTest(n = n, d = d, alpha = x, sided = sided) - power, 1e-10, 1 - 1e-10)
end

function NormTest(;
    n::Real = 0,
    d::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two")

    if sum([x == 0 for x in (n,d,alpha,power)]) != 1
        error("exactly one of n, d, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerNormTest(n = n, d = d, alpha = alpha, sided = sided)
    elseif alpha == 0.0
        alpha = alphaNormTest(n = n, d = d, power = power, sided = sided)
    elseif d == 0.0
        d = effectsizeNormTest(n = n, alpha = alpha, power = power, sided = sided)
    elseif n == 0
        n = samplesizeNormTest(d = d, alpha = alpha, power = power, sided = sided)
    end

    alt = Dict("two" => "two-sided", "less" => "less", "greater" => "greater")

    return htest(
        string("Mean power calculation for normal distribution with known variance"),
        OrderedDict(
            "n" => n,
            "d" => d,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[sided],
            "note" => "")
        )
end

#@time print(NormTest(n=0,d=.2))

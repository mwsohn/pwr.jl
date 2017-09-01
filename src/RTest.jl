# power analysis functions for mean for normal distribution with known variance

function powerRTest(;
    r::Float64 = 0.0,
    n::Real = 0,
    alpha = 0.05,
    sided::String = "two")

    check_args(r=r,n=n,alpha=alpha)

    if sided == "less"
        tside = 1
        r = -abs(r)
    elseif sided == "two"
        tside = 2
        r = abs(r)
    elseif sided == "greater"
        tside = 3
    end

    ttt = cquantile(TDist(n-2),alpha)
    rc = sqrt(ttt^2 / (ttt^2 + n - 2))
    zr = atanh(r) + r / (2 * (n - 1))
    zrc = atanh(rc)

    if tside == 1
        return cdf(Normal(),(zr - zrc) * sqrt(n-3))
    elseif tside == 2
        return cdf(Normal(),(zr - zrc) * sqrt(n-3))
            + cdf(Normal(),(-zr - zrc) * sqrt(n-3))
    elseif tside == 3
        return cdf(Normal(),(zr - zrc) * sqrt(n-3))
    end
end

function samplesizeRTest(;
    r::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two")

    check_args(r=r,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerRTest(n = x, r = r, alpha = alpha, sided = sided) - power, 3.0, 10.0^7))
end

function effectsizeRTest(;
    n::Int64 = 0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x -> powerRTest(n = n, r = x, alpha = alpha, sided = sided) - power,.001,100)
end

function alphaRTest(;
    n = 0,
    r = 0.0,
    power = 0.8,
    sided::String = "two"
    )

    check_args(r=r,n=n,power=power)

    return fzero(x->powerRTest(n = n, r = r, alpha = x, sided = sided) - power, 1e-10, 1 - 1e-10)
end

function RTest(;
    n::Real = 0,
    r::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8,
    sided::String = "two")

    if sum([x == 0 for x in (n,r,alpha,power)]) != 1
        error("exactly one of n, r, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerRTest(n = n, r = r, alpha = alpha, sided = sided)
    elseif alpha == 0.0
        alpha = alphaRTest(n = n, r = r, power = power, sided = sided)
    elseif r == 0.0
        r = effectsizeRTest(n = n, alpha = alpha, power = power, sided = sided)
    elseif n == 0
        n = samplesizeRTest(r = r, alpha = alpha, power = power, sided = sided)
    end

    alt = Dict("two" => "two-sided", "less" => "less", "greater" => "greater")

    return htest(
        string("Approximate correlation power calculation (arctangh transformation)"),
        OrderedDict(
            "n" => n,
            "r" => r,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[sided])
        )
end
# 
# @time print(RTest(n=193,r=.2,power=0.0))
# @time print(RTest(n=193,r=.2,power=0.0,sided="less"))
# @time print(RTest(n=193,r=.2,power=0.0,sided="greater"))

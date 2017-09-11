# power analysis functions for proportions

using Distributions, Roots

function powerPTest(;
    h::Real = 0,
    n::Real = 0,
    alpha::Float64 = 0.05,
    alternative::String = "two.sided"
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
        error("`alternative` must be `less`, `two.sided`, or `greater`.")
    end

    if tside == 1
        return cdf(Normal(),quantile(Normal(),alpha) - h*sqrt(n))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - h*sqrt(n)) + cdf(Normal(),quantile(Normal(),alpha/2) - h*sqrt(n))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha) - h*sqrt(n))
    end
    error("internal error")
end

#powerPTest(h=.3,n=100)

function samplesizePTest(;
    h::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(h=h,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerPTest(h = h, n = x, alpha = alpha, alternative = alternative) - power, 2 + 1e-10, 1e+09))
end

# samplesizePTest(h = .3) # 175

function effectsizePTest(;
    n::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(n=n,alpha=alpha,power=power)

    return fzero(x->powerPTest(h = x, n = n, alpha = alpha, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

# effectsizePTest(n=175) # .3

function alphaPTest(;
    h::Real = 0,
    n::Real = 0,
    power::Float64 = 0.8,
    alternative::String = "two"
    )

    check_args(h=h,n=n,power=power)

    return fzero(x->powerPTest(h = h, n = n, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

# alphaPTest(h=.3,n=175) # 0.0495

function PTest(;
    h::Real = 0,
    n::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    alternative::String = "two.sided"
    )
    if sum([x == 0 for x in (h,n,alpha,power)]) != 1
        error("exactly one of `h`, `n`, `power`, and `alpha` must be zero")
    end

    if power == 0.0
        power = powerPTest(h = h, n = n, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaPTest(h = h, n = n, power = power, alternative = alternative)
    elseif h == 0
        h = effectsizePTest(n = n, alpha = alpha, power = power, alternative = alternative)
    elseif n == 0
        n = samplesizePTest(h = h, alpha = alpha, power = power, alternative = alternative)
    end

    alt = Dict("two" => "two-sided", "two.sided" => "two-sided","less" => "less", "greater" => "greater")

    return htest(
        string("Proportion power calculation for binomial distribution (arcsine transformation)"),
        OrderedDict(
            "h" => h,
            "n" => n,
            "alpha" => alpha,
            "power" => power,
            "alternative" => alt[alternative]
            )
        )
end

#PTest(h=.0,n=100,power = 0.8)

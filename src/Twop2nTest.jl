# power analysis functions for difference in two proportions in two different sample sizes

function powerTwop2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Real = 0.05,
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

function samplesizeTwop2nTest(;
    h::Real = 0,
    n1::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(h=h,n1=n1,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerTwop2nTest(h = h, n1 = n1, n2 = x, alpha = alpha, alternative = alternative) - power, 2 + 1e-10, 1e+09))
end

function effectsizeTwop2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(n1=n1,n2=n2,alpha=alpha,power=power)

    return fzero(x->powerTwop2nTest(h = x, n1 = n1, n2 = n2, alpha = alpha, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

function alphaTwop2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    power::Real = 0.0,
    alternative::String = "two"
    )

    check_args(h=h,n1=n1,n2=n2,power=power)

    return fzero(x->powerTwop2nTest(h = h, n1 = n1, n2 = n2, alpha = x, alternative = alternative) - power, 1e-10, 1 - 1e-10)
end

function Twop2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Real = 0.05,
    power::Real = 0.0,
    alternative::String = "two"
    )
    if sum([x == 0 for x in (h,n1,n2,alpha,power)]) != 1
        error("exactly one of `h`, `n2`, `power`, and `alpha` must be zero")
    end

    if n1 == 0
        error("`n1` cannot be zero")
    end

    if power == 0.0
        power = powerTwop2nTest(h = h, n1 = n1, n2 = n2, alpha = alpha, alternative = alternative)
    elseif alpha == 0.0
        alpha = alphaTwop2nTest(h = h, n1 = n1, n2 = n2, power = power, alternative = alternative)
    elseif h == 0
        h = effectsizeTwop2nTest(n1 = n1, n2 = n2, alpha = alpha, power = power, alternative = alternative)
    elseif n2 == 0
        n2 = samplesizeTwop2nTest(h = h, n1 = n1, alpha = alpha, power = power, alternative = alternative)
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

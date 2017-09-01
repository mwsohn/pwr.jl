# power analysis functions for difference in two proportions in two different sample sizes

function power2p2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    sided::String = "two"
    )

    check_args(h=h,n1=n1,n2=n2,alpha=alpha)

    if sided == "less"
        tside = 1
    elseif sided == "two"
        tside = 2
        h = abs(h)
    elseif sided == "greater"
        tside = 3
    else
        error("`sided` must be `less`, `two`, or `greater`.")
    end

    if tside == 1
        return cdf(Normal(),quantile(Normal(),alpha) - h*sqrt(n1*n2/(n1 + n2)))
    elseif tside == 2
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - h*sqrt(n1*n2/(n1 + n2))) + cdf(Normal(),quantile(Normal(),alpha) - h*sqrt(n1*n2/(n1 + n2)))
    elseif tside == 3
        return ccdf(Normal(),cquantile(Normal(),alpha/2) - h*sqrt(n1*n2/(n1 + n2)))
    end
    error("internal error")
end

function samplesize2p2nTest(;
    h::Real = 0,
    n1::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    sided::String = "two"
    )

    check_args(h=h,n1=n1,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->power2p2nTest(h = h, n1 = n1, n2 = x, alpha = alpha, sided = sided) - power, 2 + 1e-10, 1e+09))
end

#samplesize2p2nTest(h=.3,n1=100)

function effectsize2p2nTest(;
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    sided::String = "two"
    )

    check_args(n1=n1,n2=n2,alpha=alpha,power=power)

    return fzero(x->power2p2nTest(h = x, n1 = n1, n2 = n2, alpha = alpha, sided = sided) - power, 1e-10, 1 - 1e-10)
end

# effectsizeTwopTest(n=175) # .3

function alpha2p2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    power::Float64 = 0.8,
    sided::String = "two"
    )

    check_args(h=h,n1=n1,n2=n2,power=power)

    return fzero(x->power2p2nTest(h = h, n1 = n1, n2 = n2, alpha = x, sided = sided) - power, 1e-10, 1 - 1e-10)
end

# alphaTwopTest(h=.3,n=175) # 0.0495

function Twop2nTest(;
    h::Real = 0,
    n1::Real = 0,
    n2::Real = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8,
    sided::String = "two"
    )
    if sum([x == 0 for x in (h,n1,n2,alpha,power)]) != 1
        error("exactly one of `h`, `n1`, `n2`, `power`, and `alpha` must be zero")
    end

    if n1 == 0
        n1 = n2
        n2 = 0
    end

    if power == 0.0
        power = power2p2nTest(h = h, n1 = n1, n2 = n2, alpha = alpha, sided = sided)
    elseif alpha == 0.0
        alpha = alpha2p2nTest(h = h, n1 = n1, n2 = n2, power = power, sided = sided)
    elseif h == 0
        h = effectsize2p2nTest(n1 = n1, n2 = n2, alpha = alpha, power = power, sided = sided)
    elseif n2 == 0
        n = samplesize2p2nTest(h = h, n1 = n1, alpha = alpha, power = power, sided = sided)
    end

    note = "different sample sizes"

    return htest(
        "Difference of proportion power calculation for binomial distribution (arcsine transformation)",
        OrderedDict(
            "h" => h,
            "n1" => n1,
            "n2" => n2,
            "alpha" => alpha,
            "power" => power,
            "note" => note)
        )
end

#Twop2nTest(h=.3,n1=100,n2=100,power = 0.0)

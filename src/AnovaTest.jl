# power analysis functions for ANOVA

# k = number of groups
# n = number of observations per group
# f = effect size
# alpha = Type I error
# power = 1 - Type II error
function powerAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f::Float64 = 0.0,
    alpha = 0.05)

    check_args(k=k,n=n,f=f,alpha=alpha)

    λ = k*n*f^2
    return ccdf(NoncentralF(k-1,(n-1)*k,λ),cquantile(FDist(k-1,(n-1)*k),alpha))
end

#powerAnovaTest(k=2,n=100,f=.3)

function samplesizeAnovaTest(;
    k::Real = 0,
    f::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8)

    check_args(k=k,f=f,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerAnovaTest(k = k, n = x, f = f, alpha = alpha) - power, 2.0, 10.0^7))
end

function groupsizeAnovaTest(;
    n::Real = 0,
    f::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8)

    check_args(n=n,f=f,alpha=alpha,power=power)

    return ceil(Int64,fzero(x->powerAnovaTest(k = k, n = x, f = f, alpha = alpha) - power, 2.0, 10.0^7))
end

function effectsizeAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    alpha = 0.05,
    power = 0.8)

    check_args(k=k,n=n,alpha=alpha,power=power)

    return fzero(x->powerAnovaTest(k = k, n = n, f = x, alpha = alpha) - power, 1e-07, 1e+07)
end

function alphaAnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f = 0.0,
    power = 0.8)

    check_args(k=k,f=f,n=n,power=power)

    return fzero(x->powerTTest(n = n, d = d, f = f, alpha = x) - power, 1e-10, 1 - 1e-10)
end

function AnovaTest(;
    k::Real = 0,
    n::Real = 0,
    f::Float64 = 0.0,
    alpha = 0.05,
    power = 0.8)

    if sum([x == 0 for x in (k,n,f,alpha,power)]) != 1
        error("exactly one of k, n, f, power, and alpha must be zero")
    end

    if power == 0.0
        power = powerAnovaTest(k = k, n = n, f = f, alpha = alpha)
    elseif alpha == 0.0
        alpha = alphaAnovaTest(k = k, n = n, f = f, power = power)
    elseif f == 0.0
        f = effectsizeAnovaTest(k = k, n = n, alpha = alpha, power = power)
    elseif n == 0
        n = samplesizeAnovaTest(k = k, f = f, alpha = alpha, power = power)
    elseif k == 0
        k = groupsizeAnovaTest(n = n, f = f, alpha = alpha, power = power)
    end

    alt = Dict("two" => "two-sided", "less" => "less", "greater" => "greater")

    note = "`n` is the number in each group"

    return htest(
        string("Balanced one-way analysis of variance power calculation"),
        OrderedDict(
            "k" => k,
            "n" => n,
            "f" => f,
            "alpha" => alpha,
            "power" => power,
            "note" => note)
        )
end


#print(AnovaTest(k=2,n=100,f=.3,power=0.0))

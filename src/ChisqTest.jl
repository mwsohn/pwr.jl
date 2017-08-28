# power analysis functions for chi square test

using Distributions, Roots

function powerChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    alpha::Float64 = 0.05
    )

    if w <= 0
        error("w must be positive")
    end

    if N < 1
        error("Number of observations `N` must be at least 1")
    end

    if alpha == 0.0
        error("`alpha` must be a number in [0,1]")
    end

    if df == 0
        error("`df` must be at least 1")
    end

    λ = N * w^2
    return ccdf(NoncentralChisq(df,λ), cquantile(Chisq(df),alpha))
end

function samplesizeChisqTest(;
    w::Real = 0,
    df::Int64 = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8
    )
    return fzero(x->powerChisqTest(w = w, N = x, df = df, alpha = alpha) - power, 2.0, 10.0^7)
end

function effectsizeChisqTest(;
    N::Real = 0,
    df::Int64 = 0,
    alpha::Float64 = 0.05,
    power::Float64 = 0.8
    )
    return fzero(x->powerChisqTest(w = x, N = N, df = df, alpha = alpha) - power, 1e-10, 1e+9)
end

function alphaChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    power::Float64 = 0.8
    )
    return fzero(x->powerChisqTest(w = w, N = N, df = df, alpha = x) - power, 1e-10, 1-1e-10)
end

function ChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    alpha = 0.05,
    power::Float64 = 0.8
    )
    if sum([x == 0 for x in (w,N,alpha,power)]) != 1
        error("exactly one of `w`, `N`, `power`, and `alpha` must be zero")
    end

    if power == 0.0
        power = powerChisqTest(w = w, N = N, df = df, alpha = alpha)
    elseif alpha == 0.0
        alpha = alphaChisqTest(w = w, N = N, df = df, power = power)
    elseif w == 0
        w = effectsizeChisqTest(N = N, df = df, alpha = alpha, power = power)
    elseif N == 0
        N = samplesizeChisqTest(w = w, df = df, alpha = alpha, power = power)
    end

    println("\nChi-square test power calculation\n")
    @printf("%13s = %.6f\n","w",w)
    @printf("%13s = %d\n","N",N)
    @printf("%13s = %d\n","df",df)
    @printf("%13s = %.6f\n","alpha",alpha)
    @printf("%13s = %.6f\n","power",power)
    println("\nNOTE: N is the number of observations")
end

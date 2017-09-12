# power analysis functions for chi square test

function powerChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05
    )

    check_args(w = w, N = N, df = df, alpha = alpha)

    λ = N * w^2
    return ccdf(NoncentralChisq(df,λ), cquantile(Chisq(df),alpha))
end

function samplesizeChisqTest(;
    w::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05,
    power::Real = 0.0
    )

    check_args(w = w, df = df, alpha = alpha, power = power)

    return fzero(x->powerChisqTest(w = w, N = x, df = df, alpha = alpha) - power, 2.0, 10.0^7)
end

function effectsizeChisqTest(;
    N::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05,
    power::Real = 0.0
    )

    check_args(N = N, df = df, alpha = alpha, power = power)

    return fzero(x->powerChisqTest(w = x, N = N, df = df, alpha = alpha) - power, 1e-10, 1e+9)
end

function alphaChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    power::Real = 0.0
    )

    check_args(w = w, N = N, df = df, power = power)

    return fzero(x->powerChisqTest(w = w, N = N, df = df, alpha = x) - power, 1e-10, 1-1e-10)
end

function ChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05,
    power::Real = 0.0
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

    note = "`N` is the number of observations"

    return htest(
        string("Chi-square test power calculation"),
        OrderedDict(
            "w" => w,
            "N" => N,
            "df" => df,
            "alpha" => alpha,
            "power" => power,
            "note" => note)
        )

end

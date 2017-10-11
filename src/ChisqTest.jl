# power analysis functions for chi square test

"""
    powerChisqTest(w::Real = 0, N::Real = 0, df::Int64 = 0.0, alpha::Real = 0.05)

Estimates power of a sample with `N` observations and `df` degrees of freedom
to detect an effect size `w` at type I error = `alpha` (default: 0.05).
The effect size `w` is the square root of the sum of the squared differences between pᵢ and p₀
divided by p₀ in each cell, where pᵢ is the proportion under the alternative hypothesis
and p₀ is the proportion under the null hypothesis. Two functions `ESw1(P₀,P₁)` and `ESw2(P)`
are available for you to estimate `w`.
"""
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

"""
    samplesizeChisqTest(w::Real = 0, df::Int64 = 0.0, alpha::Real = 0.05, power::Float64 = 0.8)

Estimates the minimum sample size required to detect an effect size `w` in a table
with `df` degrees of freedom at power > `power` (default: 0.8) and type I error = `alpha` (default: 0.05).
The effect size `w` is the square root of the sum of the squared differences between pᵢ and p₀
divided by p₀ in each cell, where pᵢ is the proportion under the alternative hypothesis
and p₀ is the proportion under the null hypothesis. Two functions `ESw1(P₀,P₁)` and `ESw2(P)`
are available for you to estimate `w`.
"""
function samplesizeChisqTest(;
    w::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05,
    power::Real = 0.8
    )

    check_args(w = w, df = df, alpha = alpha, power = power)

    return ceil(Int64,fzero(x->powerChisqTest(w = w, N = x, df = df, alpha = alpha) - power, 2.0, 10.0^7))
end

"""
    effectsizeChisqTest(N::Real = 0, df::Int64 = 0.0, alpha::Real = 0.05, power::Float64 = 0.8)

Estimates effect size the sample with `N` observations with `df` degrees of freedom
can detect with power > `power` (default: 0.8) at type I error = `alpha` (default: 0.05).
"""
function effectsizeChisqTest(;
    N::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.05,
    power::Real = 0.0
    )

    check_args(N = N, df = df, alpha = alpha, power = power)

    return fzero(x->powerChisqTest(w = x, N = N, df = df, alpha = alpha) - power, 1e-10, 1e+9)
end

"""
    alphaChisqTest(w::Real = 0, N::Real = 0, df::Int64 = 0.0, power::Float64 = 0.8)

Estimates the probability of type I error (alpha) that the sample `N` has in detecting
an effect size `w` in a table with `df` degrees of freedom at power > `power` (default: 0.8).
The effect size `w` is the square root of the sum of the squared differences between pᵢ and p₀
divided by p₀ in each cell, where pᵢ is the proportion under the alternative hypothesis
and p₀ is the proportion under the null hypothesis. Two functions `ESw1(P₀,P₁)` and `ESw2(P)`
are available for you to estimate `w`.
"""
function alphaChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    power::Real = 0.0
    )

    check_args(w = w, N = N, df = df, power = power)

    return fzero(x->powerChisqTest(w = w, N = N, df = df, alpha = x) - power, 1e-10, 1-1e-10)
end

"""
    pwr.ChisqTest(w::Real = 0, N::Real = 0, df::Int64 = 0.0, alpha::Float64 = 0.0, power::Float64 = 0.0)

Estimates one of the test parameters such as sample size (`N`),
effect size (`w`), degree of freedom (`df`), type I error (`alpha`),
and power (`power`). The parameter to be estimated must be set to zero (default).
The effect size `w` is the square root of the sum of the squared differences between pᵢ and p₀
divided by p₀ in each cell, where pᵢ is the proportion under the alternative hypothesis
and p₀ is the proportion under the null hypothesis. Two functions `ESw1(P₀,P₁)` and `ESw2(P)`
are available for you to estimate `w`.
"""
function ChisqTest(;
    w::Real = 0,
    N::Real = 0,
    df::Int64 = 0,
    alpha::Real = 0.0,
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

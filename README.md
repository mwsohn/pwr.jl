# pwr
Julia port of `pwr` package in R. This package was originally created by Stephane Champely, from the University of Lyon.

## Installation
To install `pwr`, use the following:

```jldoctest
Pkg.clone("https://github.com/mwsohn/pwr.jl")
```

To use the package, start by `using pwr` in your session.

## Tests for power calculation

All functions in the original R package was ported. For each family for tests,
four functions were made available to compute `power`, `samplesize`, `effectsize`,
and `alpah (Type I error)`. These functions were all exported.

### ANOVA

Options:

- `k`: number of groups
- `n` = number of observations per group
- `f` = effect size
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)

Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For example, for `powerAnovaTest()`, do not specify `power` option: powerAnovaTest(k=2,n=200,f=.3) will
compute power for a sample of 200 in each group with an effect size of .3. For `pwr.AnovTest()`, exactly
one of the options needs to be set to zero.

- `powerAnovaTest()`
- `samplesizeAnovaTest()`
- `effectsizeAnovaTest()`
- `alphaAnovaTest()`
- `pwr.AnovaTest()`

#### Examples

```jldoctest
julia> p = powerAnovaTest(n=100,k=2,f=.2)
0.8036475048589252

julia> tst = pwr.AnovaTest(n=100,k=2,f=.2,power=0.0)

Balanced one-way analysis of variance power calculation

            k = 2
            n = 100
            f = 0.2
        alpha = 0.05
        power = 0.8036475048589252

NOTE: `n` is the number in each group


julia> plot(tst)
```

### T-test

Options:

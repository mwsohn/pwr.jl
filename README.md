# pwr
Julia port of `pwr` package in R. This package was originally created by Stephane Champely, from the University of Lyon.

## Installation
To install `pwr`, use the following:

```jldoctest
Pkg.clone("https://github.com/mwsohn/pwr.jl")
```

To use the package, start by `using pwr` in your session.

## Tests for power calculation

All functions in the original R package was ported. For each family of tests,
four functions additional are available that start with `power`, `samplesize`, `effectsize`,
and `alpha` (Type I error). These functions were all exported.

### pwr.AnovaTest

Power calculations for balanced one-way analysis of variance tests

#### Options:

- `k`: number of groups
- `n` = number of observations per group
- `f` = effect size
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error

#### Functions:

All options listed above must be specified as keyword argments except for the parameter to be estimated.
For example, for `powerAnovaTest()`, do not specify `power` option: powerAnovaTest(k=2,n=200,f=.3) will
compute power for a sample of 200 in each group with an effect size of .3. For `pwr.AnovaTest()`, exactly
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

### pwr.TTest

Power calculations for t-test of means (one sample, two samples and paired samples)

#### Options:

- `d`: effect size
- `n` = number of observations per group
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error
- `sampletype` = "onesample" for one sample t-test or "twosample" for two sample t-test (default: "onesample")
- `alternative` = "less" for testing μ₁ < μ₂, "two-sided" for μ₁ = μ₂, and "greater" for μ₁ > μ₂ (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For example, for `powerTTest()`, do not specify `power` option: powerTTest(k=2,n=200,d=.3) will
compute power for a sample of 200 in each group with an effect size of .3. For `pwr.TTest()`, exactly
one of the options `d`, `n`, `alpha`, and `power` needs to be set to zero.

- `powerTTest()`
- `samplesizeTTest()`
- `effectsizeTTest()`
- `alphaTTest()`
- `pwr.TTest()`

#### Examples

```jldoctest
julia> powerTTest(n=100,d=.3, alternative="two-sided")
0.8439471027376636

julia> samplesizeTTest(d=.3, power = 0.8, alternative="two-sided")
90

julia> tst = pwr.TTest(n=100,d=.3,power=0.0, sampletype = "onesample", alternative="two-sided")
One-sample t-test power calculation

            n = 100
            d = 0.3
        alpha = 0.05
        power = 0.8439471027376636
   sampletype = One-sample
  alternative = two-sided

NOTE: `n` is number in each group

julia> plot(tst)
```

### T2nTest

Compute power of tests or determine parameters to obtain target power (similar to TTest)

#### Options

- `n1` = Number of observations in the first sample
- `n2` = Number of observations in the second sample
- `d` = Effect size
- `alpha` = Type I error
- `power` = Power of test (1 minus Type II error)
- `alternative` = "less" for Hₐ:μ₁ < μ₂, "two-sided" for Hₐ:μ₁ = μ₂, and "greater" for Hₐ:μ₁ > μ₂ (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For example, for `powerTTest()`, do not specify `power` option: powerTTest(k=2,n=200,d=.3) will
compute power for a sample of 200 in each group with an effect size of .3. For `pwr.TTest()`, exactly
one of the options `d`, `n`, `alpha`, and `power` needs to be set to zero.

- `powerTTest()`
- `samplesizeTTest()`
- `effectsizeTTest()`
- `alphaTTest()`
- `pwr.TTest()`

#### Examples

```jldoctest
julia> powerT2nTest(n1=100,n2=150,d=.3, alternative="two-sided")
0.6386472954918861

julia> effectsizeT2nTest(n1=100,n2=150,power = 0.8,alternative="two-sided")
0.3630907608326113

julia> pwr.T2nTest(n1=100,n2 = 0,d=.3,power=0.0, alternative="two-sided")
T-test power calculation

           n1 = 100
           n2 = 695
            d = 0.3
        alpha = 0.05
        power = 0.8
  alternative = two-sided
```

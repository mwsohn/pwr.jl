# pwr
Julia port of `pwr` package in R. This package was originally created by Stephane Champely, from the University of Lyon. This package implemented power calculations along the lines of Cohen (1988) using in particular the same notations for effect sizes. Some examples from the book are used below.

## Installation
To install `pwr`, use the following:

```jldoctest
Pkg.clone("https://github.com/mwsohn/pwr.jl")
```

To use the package, start by `using pwr` in your session.

## Tests for power calculation

This package contains functions for basic power calculations using effect sizes and notations from Cohen (1988):

- [pwr.AnovaTest](#pwranovatest): test for one-way balanced anova (ES=f)
- [pwr.ChisqTest](#pwrchisqtest): chi-squared test (ES=w)
- [pwr.F2Test](#pwrf2test): test for the general linear model (ES=f2)
- [pwr.PTest](#pwrptest): test for one proportion (ES=h)
- [pwr.RTest](#pwrrtest): correlation test (ES=r)
- [pwr.TwopTest](#pwrtwoptest): test for two proportions (ES=h)
- [pwr.Twop2nTest](#pwrtwo2ntest): test for two proportions (ES=h, unequal sample sizes)
- [pwr.TTest](#pwrttest): one sample and two samples (equal sizes) t tests for means (ES=d)
- [pwr.T2nTest](#pwrt2ntest): two samples (different sizes) t test for means (ES=d)

The following utility functions are also available:

- [cohenES](#cohenes): computing effect sizes for all the previous tests corresponding to conventional effect sizes (small, medium, large)
- [ESh](#esh): computing effect size h for proportions tests
- [ESw1](#esw1): computing effect size w for the goodness of fit chi-squared test
- [ESw2](#esw2): computing effect size w for the association chi-squared test
- [plot](#plot): plot sample sizes and their estimated power

All functions in the original R package were ported. For each family of tests,
four functions are additionally available that start with `power`, `samplesize`, `effectsize`, and `alpha` to specifically compute power, sample size, effect size, and alpha (Type I error), respectively. These functions were all exported, but the main tests were not and should be used with `pwr.` in front.

`fzero` function in the `Roots.jl` module is used to solve power equations for unknown values, so you may see errors from it. These errors may arise because the unknown value could not be computed with the given values, notably sample sizes. If you have such an error, consider increasing the sample size or the effect size.

### pwr.AnovaTest

Compute power for balanced one-way analysis of variance tests

#### Options:

- `k`: number of groups
- `n` = number of observations per group
- `f` = effect size
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)

#### Functions:

All options listed above must be specified as keyword argments
except for the parameter to be estimated. For example,
for `powerAnovaTest()`, do not specify `power` option:
powerAnovaTest(k=2,n=200,f=.3) will compute power for a sample of 200
in each group with an effect size of .3. For `pwr.AnovaTest()`, exactly
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
![Image](plots/anova.png?raw=true)

### pwr.ChisqTest

Compute power of test or determine parameters to obtain target power (same as power.anova.test)

#### Options:

- `N` = total number of observations
- `w` = effect size
- `df` = degree of freedom (depends of the chosen test)
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For example, for `powerChisqTest()`, do not specify `power` option: powerChisqTest(df=2,N=200,w=.3) will
compute power for a sample of 200 in a table with a degree of freedom 2 with an effect size of .3. For `pwr.ChisqTest()`, exactly one of the options `w`, `N`, `alpha`, and `power` needs to be set to zero.

- `powerChisqTest()`
- `samplesizeChisqTest()`
- `effectsizeChisqTest()`
- `alphaChisqTest()`
- `pwr.ChisqTest()`

#### Examples

```jldoctest
julia> powerChisqTest(w.=0.289,df=(4-1)*(3-1),N=100,alpha=0.05)
0.5518252596952123

julia> alphaChisqTest(w=0.346,df=(2-1)*(3-1),N=140,power=0.8)
0.00321629291888053

julia> tst = pwr.ChisqTest(w=0.1,df=(5-1)*(6-1),power=0.80,alpha=0.05)
Chi-square test power calculation

            w = 0.1
            N = 2097
           df = 20
        alpha = 0.05
        power = 0.8

NOTE: `N` is the number of observations

julia> plot(tst)
```
![Image](plots/chisq.png?raw=true)

### pwr.F2Test

Compute power of test or determine parameters to obtain target power.
This function provides power calculations for the general linear model.

#### Options:

- `u` = degree of freedom for numerator
- `v` = degree of freedom for denominator
- `f2` = effect size
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For `pwr.F2Test()`, exactly one of the options `u`, `v`, `f2`, `alpha`, and `power` needs to be set to zero. Plot for F2Test is not supported.

- `powerF2Test()`
- `samplesizeF2Test()`
- `effectsizeF2Test()`
- `alphaF2Test()`
- `pwr.F2Test()`

#### Examples

```jldoctest
julia> pwr.F2Test(u=5,v=89,f2=0.1/(1-0.1),alpha=0.05)
Multiple regression power calculation

            u = 5
            v = 89
           f2 = 0.11111111111111112
        alpha = 0.05
        power = 0.6735857709143758
```

### pwr.NormTest

Compute power for the mean of a normal distribution with known variance

#### Options:

- `d` = effect size (μ - μ₀)
- `n` = number of observations
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)
- `alternative` = "less" for testing μ < μ₀, "two-sided" for μ = μ₀, and "greater" for μ > μ₀ (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For `pwr.NormTest()`, exactly one of the options `d`, `n`, `alpha`, and `power` needs to be set to zero.

- `powerNormTest()`
- `samplesizeNormTest()`
- `effectsizeNormTest()`
- `alphaNormTest()`
- `pwr.NormTest()`

#### Examples

Power at μ = 105 for H₀:μ=100 vs. Hₐ:μ>100 (σ = 15) in a sample of 20 observations
at alpha = 0.05 can be caculated as follows:

```jldoctest
julia> σ = 15.
15.0

julia> c = 100
100

julia> μ = 105.
105.0

julia> d = (μ -c)/σ
0.3333333333333333

julia> t = pwr.NormTest(d=d,n=20,alpha=0.05,alternative="greater")
Mean power calculation for normal distribution with known variance

            n = 20
            d = 0.3333333333333333
        alpha = 0.05
        power = 0.4387490275410111
  alternative = greater

julia> samplesizeNormTest(d=d,power=0.8,alpha=0.05,alternative="greater")
56

julia> μ = collect(linspace(95,125,100))
100-element Array{Float64,1}:
  95.0   
  95.303
  95.6061
  95.9091
  96.2121
  96.5152
  96.8182
  97.1212
  97.4242
  97.7273
   ⋮     
 122.576
 122.879
 123.182
 123.485
 123.788
 124.091
 124.394
 124.697
 125.0   

julia> d = (μ-c)/σ
100-element Array{Float64,1}:
 -0.333333
 -0.313131
 -0.292929
 -0.272727
 -0.252525
 -0.232323
 -0.212121
 -0.191919
 -0.171717
 -0.151515
  ⋮       
  1.50505
  1.52525
  1.54545
  1.56566
  1.58586
  1.60606
  1.62626
  1.64646
  1.66667

julia> power = [powerNormTest(d=x,n=20,alternative="greater") for x in d]
100-element Array{Float64,1}:
 0.000857615
 0.00116255
 0.00156399
 0.00208816
 0.00276704
 0.00363915
 0.00475039
 0.0061548  
 0.00791534
 0.0101044  
 ⋮          
 1.0        
 1.0        
 1.0        
 1.0        
 1.0        
 1.0        
 1.0        
 1.0        
 1.0        
julia> plot(d,power,ylim=[0.,1.],legend=false,ylabel="Test Power = 1 - \\beta", xlabel = "Effect Size")

julia> hline!([0.05,0.80])
```
![Image](plots/f2_1.png?raw=true)

```jldoctest
julia> plot(d,[powerNormTest(d=x,n=20,alpha=0.05,alternative="two.sided") for x in d],ylim=[0,1],legend=false)

julia> hline!([0.05,0.8])
```
![Image](plots/f2_2.png?raw=true)

### pwr.PTest

Compute power for proportion tests (one sample). These calculations use arcsine transformation of the proportion (see Cohen (1988)). Use `ESh()` (see below) to compute the effect size from two proportions.

#### Options:

- `h` = effect size
- `n` = number of observations per group
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)
- `alternative` = "less" for testing p₁ < p₂, "two-sided" for p₁ = p₂, and "greater" for p₁ > p₂ (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For `pwr.PTest()`, exactly one of the options `h`, `n`, `alpha`, and `power` needs to be set to zero.

- `powerPTest()`
- `samplesizePTest()`
- `effectsizePTest()`
- `alphaPTest()`
- `pwr.PTest()`

#### Examples

```jldoctest
julia> h = ESh(0.5,0.4)
0.20135792079033088

julia> powerPTest(h=h,n=60,alpha=0.05,alternative="two.sided")
0.3447014091272134

julia> tst = pwr.PTest(h=0.2,power=0.95,alpha=0.05,alternative="two.sided")
Proportion power calculation for binomial distribution (arcsine transformation)

            h = 0.2
            n = 325
        alpha = 0.05
        power = 0.95
  alternative = two-sided

julia> plot(tst)
```
![Image](plots/p.png?raw=true)


### pwr.RTest

Compute power for correlation test. These calculations use the Z’ transformation of correlation coefficient : Z’=arctanh(r)+r/(2*(n-1)) (see Cohen (1988) p.546).

#### Options:

- `n` = number of observations per group
- `r`= linear correlation coefficient
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)
- `alternative` = "less", "two-sided", or "greater" (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For `pwr.RTest()`, exactly one of the options `n`, `r`, `alpha`, and `power` needs to be set to zero.

- `powerRTest()`
- `samplesizeRTest()`
- `effectsizeRTest()`
- `alphaRTest()`
- `pwr.RTest()`

#### Examples

```jldoctest
julia> pwr.r.test(r=0.3,n=50,sig.level=0.05,alternative="two.sided")
0.571535679145053

julia> pwr.r.test(r=0.3,n=50,sig.level=0.05,alternative="greater")
0.6911394853796565

julia> samplesizeRTest(r=0.3,power=0.80,alpha=0.05,alternative="two.sided")
85

julia> samplesizeRTest(r=0.5,power=0.80,alpha=0.05,alternative="two.sided")
29

julia> samplesizeRTest(r=0.5,power=0.80,alpha=0.05,alternative="two.sided")
782

julia> tst = pwr.RTest(r=.3,power=0.8, alternative="two.sided")
Approximate correlation power calculation (arctangh transformation)

            n = 85
            r = 0.3
        alpha = 0.05
        power = 0.8
  alternative = two-sided

julia> plot(tst)
```
![Image](plots/r.png?raw=true)


### pwr.TTest

Compute power for t-test of means (one sample, two samples and paired samples)

#### Options:

- `d` = effect size
- `n` = number of observations per group
- `alpha` = Type I error (default: 0.05)
- `power` = 1 - Type II error (default: 0.8)
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

julia> tst = pwr.TTest(d=.3,power=0.8, sampletype = "onesample", alternative="two.sided")
One-sample t-test power calculation

            n = 90
            d = 0.3
        alpha = 0.05
        power = 0.8
   sampletype = One-sample
  alternative = two-sided

NOTE: `n` is number in each group

julia> plot(tst)
```
![Image](plots/t.png?raw=true)

### pwr.T2nTest

Compute power of two samples (different sizes) t-tests of means.

#### Options

- `n1` = Number of observations in the first sample
- `n2` = Number of observations in the second sample
- `d` = Effect size
- `alpha` = Type I error (default: 0.05)
- `power` = Power of test (1 minus Type II error) (default: 0.8)
- `alternative` = "less" for Hₐ:μ₁ < μ₂, "two-sided" for Hₐ:μ₁ = μ₂, and "greater" for Hₐ:μ₁ > μ₂ (default: "two-sided")

#### Functions:

All options listed above need to be specified as keyword argments except for the parameter to be estimated.
For `pwr.T2nTest()`, exactly one of the options `d`, `n2`, `alpha`, and `power` needs to be set to zero.

- `powerT2nTest()`
- `samplesizeT2nTest()`
- `effectsizeT2nTest()`
- `alphaT2nTest()`
- `pwr.T2nTest()`

#### Examples

```jldoctest
julia> powerT2nTest(n1=100,n2=150,d=.3, alternative="two-sided")
0.6386472954918861

julia> effectsizeT2nTest(n1=100,n2=150,power = 0.8,alternative="two-sided")
0.3630907608326113

julia> tst = pwr.T2nTest(n1=100,n2 = 0,d=.3,power=0.8, alternative="two.sided")
T-test power calculation

           n1 = 100
           n2 = 695
            d = 0.3
        alpha = 0.05
        power = 0.8
  alternative = two-sided

julia> plot(tst)
```
![Image](plots/t2n.png?raw=true)

## Utility functions

### cohenES

Compute effect sizes for all the previous tests corresponding to conventional effect sizes (small, medium, large)

#### Options
- `test` = `p`, `t`, `r`, `anova`, `chisq`, and `f2`
- `size` = `small`, `medium`, or `large`

#### Examples

```jldoctest
julia> cohenES(test="r",size="medium")
0.3

julia> samplesizeRTest(r=cohenES(test="r",size="medium"),power=0.8,alpha=0.05,alternative="two.sided")
85

julia> pwr.RTest(r=cohenES(test="r",size="medium"),n=0,power=0.8,alpha=0.05,alternative="two.sided")

Approximate correlation power calculation (arctangh transformation)

            n = 85
            r = 0.3
        alpha = 0.05
        power = 0.8
  alternative = two-sided

```

### ESh

Compute effect size `h` for two proportions

#### Usage

ESh(p1, p2)

#### Arguments

- `p1`
            first proportion
- `p2`
            second proportion

#### Examples

```jldoctest
julia> h = ESh(0.5,0.4)
0.20135792079033088

julia> pwr.PTest(h = h,n=60,alpha=0.05,alternative="two.sided")

Proportion power calculation for binomial distribution (arcsine transformation)

            h = 0.20135792079033088
            n = 60
        alpha = 0.05
        power = 0.3447014091272134
  alternative = two-sided
```

### ESw1

Compute effect size w for two sets of k probabilities P0 (null hypothesis) and P1 (alternative hypothesis)

#### Usage

ESw1(P0, P1)

#### Arguments

- `P0`
            A vector of the first set of k probabilities (null hypothesis)
- `P1`
            A vector of the second set of k probabilities (alternative hypothesis)

#### Examples

```jldoctest
julia> P0 = fill(0.25,4)
4-element Array{Float64,1}:
 0.25
 0.25
 0.25
 0.25

julia> P1 = vcat(0.375,fill((1-0.375)/3,3))
4-element Array{Float64,1}:
 0.375   
 0.208333
 0.208333
 0.208333

julia> ESw1(P0,P1)
0.2886751345948129

julia> pwr.ChisqTest(w = ESw1(P0,P1),N=100,df=(4-1))

Chi-square test power calculation

            w = 0.2886751345948129
            N = 100
           df = 3
        alpha = 0.05
        power = 0.6739833924326929

NOTE: `N` is the number of observations

```

### ESw2

Compute effect size w for a two-way probability table corresponding to the alternative hypothesis
in the chi-squared test of association in two-way contingency tables

#### Usage

ESw2(P)

#### Arguments

- `P`
            A vector of the first set of k probabilities (null hypothesis)

#### Examples

```jldoctest
julia> prob = [0.225 0.125 0.125 0.125; 0.16 0.16 0.04 0.04]
2×4 Array{Float64,2}:
 0.225  0.125  0.125  0.125
 0.16   0.16   0.04   0.04

julia> ESw2(prob)
0.2558646068639514

julia> pwr.ChisqTest(w = ESw2(prob),df=(2-1)*(4-1),N=200)

Chi-square test power calculation

            w = 0.2558646068639514
            N = 200
           df = 3
        alpha = 0.05
        power = 0.873322154622669

NOTE: `N` is the number of observations

```

### plot

Plot a diagram to illustrate the relationship of sample size and test power for a given set of parameters

#### Usage

plot(x; backend = gr)

#### Arguments

- `x`
            object of a returned `struct` usually created by one of the power calculation functions,
e.g., pwr.TTest()

- `backend` = `gr`, `plotly`, `plotlyjs`, and `pyplot` (default: gr)

#### Details

Power calculations for the following tests are supported: t-test (pwr.TTest(), pwr.T2nTest()), chi
squared test (pwr.ChisqTest()), one-way ANOVA (pwr.AnovaTest(), standard normal distribution
(pwr.NormTest()), pearson correlation (pwr.RTest()), proportions (pwr.PTest(), pwr.TwopTest(), pwr.Twop2nTest())).

`plot` is implemented using Plots.jl package. The default backend is `gr()`. When you want to modify the plot or save the plot into a graphics file, first use `Plots`.

#### Examples

```jldoctest
julis> using pwr, Plots

julia> tst = pwr.ChisqTest(w = .3,df=3,N=200)

Chi-square test power calculation

            w = 0.3
            N = 200
           df = 3
        alpha = 0.05
        power = 0.9590732637512994

NOTE: `N` is the number of observations

julia> plot(tst)

julia> hline!([0.8])

julia> png("chisq_with_hline.png")
```
![image](plots/chisq_with_hline.png?raw=true)

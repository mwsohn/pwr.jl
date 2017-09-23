using Plots, DataFrames, Roots

import RecipesBase.plot

function plot(ht::pwr.htest; backend::Function = gr)

    methodlist = (
        "One-sample t-test power calculation",
        "Two-sample t-test power calculation",
        "Paired t-test power calculation",
        "T-test power calculation",
        "Difference of proportion power calculation for binomial distribution (same sample)",
        "Difference of proportion power calculation for binomial distribution (different sample)",
        "Balanced one-way analysis of variance power calculation",
        "Chi-square test power calculation",
        "Mean power calculation for normal distribution with known variance",
        "Proportion power calculation for binomial distribution (arcsine transformation)",
        "Approximate correlation power calculation (arctangh transformation)"
    )

    if backend in (gr, pyplot, plotly, plotlyjs) == false
        error(backend," is not supported")
    end

    if backend in (gr, pyplot)
        sepstr = "\n"
    else
        sepstr = "<br>"
    end

    ytitle = "Test Power = 1 - β"
    if backend == gr
        ytitle = "Test Power = 1 - \\beta"
    end

    d = ht.d

    if !(ht.title in methodlist)
        warn("The method `",ht.title,"` is not supported. Supported methods include: \n\t",join(methodlist,"\n\t"))
        exit(0)
    end

    breaks = 20

    # case: One-sample, Two-sample or Paired t test
    if ht.title in ("One-sample t-test power calculation",
        "Two-sample t-test power calculation",
        "Paired t-test power calculation")

        n = d["n"]
        n_upper = ceil(Int64,max(n*1.5, n+30))
        n_increment = ceil(Int64,(n_upper - 10)/breaks)
        sample_sizes = collect(10:n_increment:n_upper)

        # obtain power at each sample size
        power = [ss->powerTTest(n=ss,d=d["d"],alpha = d["alpha"],sampletype=d["sampletype"],alternative=d["alternative"]) for ss in sample_sizes]

        # create labels
        title = ht.title
        legend_string = string(
                    "tails = ", d["alternative"],sepstr,
                    "effect size d = ", d["d"],sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string(
                    "optimal sample size",sepstr,
                    "n = ", ceil(Int64,n),sepstr,d["note"])

    # case: Two-sample t test with n1 and n2
    elseif ht.title == "T-test power calculation"
        n = d["n1"] + d["n2"]
        n_upper = ceil(Int64,max(n*1.5, n+30)) # upper at least 30 above n
        n_rel = d["n1"] / n # relative sample size; will be kept constant in claculations

        # generate data
        n_increment = ceil(Int64,(n_upper - 10)/breaks)
        n_lower = ceil(Int64,fzero(x->(x*n_rel) - 2.1,2,d["n2"]))
        sample_sizes = collect(n_lower:n_increment:n_upper)

        function pwrt2test(ss)
            na = ceil(Int64,ss*n_rel)
            nb = ss - na
            if (na < 2 || nb < 2)
                return 0.0
            end
            return powerT2nTest(n1=na, n2=nb, d=d["d"], alpha = d["alpha"], alternative = d["alternative"])
        end

        power = [pwrt2test(ss) for ss in sample_sizes]

        # create labels
        title = ht.title
        legend_string = string(
                    "tails = ", d["alternative"], sepstr,
                    "effect size d = ", d["d"], sepstr,
                    "alpha = ", d["alpha"], sepstr,
                    "n1/n2 = ",round(n_rel,2))
        optimal_string = string("optimal sample size",sepstr,
                    "n = ",  d["n1"], " + ", d["n2"], " = ", n)

    # case: Difference of proportion (same sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (same sample)"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ ss-> powerTwopTest(n=ss, h=d["h"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        title = string("Difference of proportion power calculation",sepstr,"for binomial distribution (same sample)")
        legend_string = string(
                    "tails = ", d["alternative"], sepstr,
                    "effect size h = ", d["h"], sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size", sepstr,
                    "n = ", ceil(Int64,n), sepstr,d["note"])

    # case: difference of proportion (different sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (different sample)"
        n = d["n1"] + d["n2"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_rel = d["n1"] / n # relative sample size; will be kept constant in claculations
        n_increment = ceil(Int64,(n_upper - 10)/breaks)
        n_lower = ceil(Int64,fzero(x->(x*n_rel) - 2.1,2,d["n2"]))

        # generate data
        sample_sizes = collect(n_lower:n_increment:n_upper)

        function pwr2p2ntest(ss)
            na = ceil(Int64,ss*n_rel)
            nb = ss - na
            if (na < 2 || nb < 2)
                return (0.0)
            end
            return powerTwop2nTest(n1=na, n2=nb, h=d["h"], alpha = d["alpha"], alternative = d["alternative"])
        end

        power = [ss->pwr2p2ntest(ss) for ss in sample_sizes]

        # create labels
        title = string("Difference of proportion power calculation",sepstr,"for binomial distribution (different sample)")
        legend_string = string(
                    "tails = ", d["alternative"],sepstr,
                    "effect size h = ", d["h"], sepstr,
                    "alpha = ", d["alpha"], sepstr,
                    "n1/n2 = ",round(n_rel,2))
        optimal_string = string("optimal sample size", sepstr,
                    "n = ",  d["n1"], " + ", d["n2"], " = ", n)

    # case: ANOVA
    elseif ht.title== "Balanced one-way analysis of variance power calculation"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerAnovaTest(n=ss, k=d["k"], f=d["f"], alpha = d["alpha"]) for ss in sample_sizes]

        # create labels
        title = ht.title
        legend_string = string(
                    "groups k = ", d["k"], sepstr,
                    "effect size f = ", d["f"], sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size",sepstr,
                    "n = ", ceil(Int64,n),sepstr,d["note"])

    # case: Chi Squared
    elseif ht.title == "Chi-square test power calculation"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerChisqTest(N=ss, w=d["w"], alpha = d["alpha"], df=d["df"]) for ss in sample_sizes]

        # create labels
        title = ht.title
        legend_string = string(
                    "df = ", d["df"], sepstr,
                    "effect size w = ", d["w"], sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size", sepstr,
                    "N = ", ceil(Int64,n),sepstr,d["note"])

    # case: Normal distribution
    elseif ht.title == "Mean power calculation for normal distribution with known variance"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerNormTest(n=ss, d=d["d"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        title = string("Mean power calculation for normal distribution",sepstr,"with known variance")
        legend_string = string(
                    "tails = ", d["alternative"],sepstr,
                    "effect size d = ", d["d"],sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size", sepstr,
                    "n = ", ceil(Int64,n),sepstr,d["note"])

    # case: proportion
    elseif ht.title == "Proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerPTest(n=ss, h=d["h"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        title = string("Proportion power calculation",sepstr,"for binomial distribution (arcsine transformation)")
        legend_string = string(
                    "tails = ", d["alternative"],sepstr,
                    "effect size h = ", d["h"],sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size",sepstr,
                    "n = ", ceil(Int64,n))

    # case: correlation
    elseif ht.title == "Approximate correlation power calculation (arctangh transformation)"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerRTest(n=ss, r=d["r"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        title = string("Approximate correlation power calculation",sepstr,"(arctangh transformation)")
        legend_string = string(
                    "tails = ", d["alternative"],sepstr,
                    "effect size r = ", d["r"],sepstr,
                    "alpha = ", d["alpha"])
        optimal_string = string("optimal sample size",sepstr,
                    "n = ", ceil(Int64,n))
    end

    # use DataFrame
    len = length(sample_sizes)
    if len < 2
        error("cannot display power for this test")
    end

    # select the backend
    backend()

    n_lower = sample_sizes[1]
    xlim_lower = sample_sizes[len] / 80

    # plot with title and x-axis and y-axis labels
    plot(sample_sizes,
        power,
        linecolor = :blue,
        marker=(2,.5,:circle,:blue),
        #left_margin = 8mm,
        top_margin = 8mm,
        right_margin = 4mm,
        title = title,
        yaxis = (ytitle,(-0.03,1.03),0.0:0.2:1.0),
        xaxis = ("Sample Size",(n_lower - xlim_lower,n_upper)),
        legend=false,
        annotations=([((n_lower+xlim_lower),0.99,text(legend_string,9,:blue,:left,:top)),
            (sample_sizes[floor(Int64,len/2)],0.06,text(optimal_string,9,:red,:left,:bottom))]))
    vline!([n],line=(1,:dot,.8,:red))

end

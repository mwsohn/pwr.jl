using pwr, DataFrames, Plots

import RecipesBase.plot

function plot(ht::pwr.htest)

    methods = (
        "One-sample t-test power calculation",
        "Two-sample t-test power calculation",
        "Paired t-test power calculation",
        "T-test power calculation",
        "Difference of proportion power calculation for binomial distribution (arcsine transformation)",
        "Difference of proportion power calculation for binomial distribution (arcsine transformation)",
        "Balanced one-way analysis of variance power calculation",
        "Chi squared power calculation",
        "Mean power calculation for normal distribution with known variance",
        "proportion power calculation for binomial distribution (arcsine transformation)",
        "Approximate correlation power calculation (arctangh transformation)"
    )

    d = ht.d

    if !(ht.title in methods)
        warn("The method `",ht.title,"` is not supported. Supported methods include: \n\t",join(methods,"\n\t"))
        exit(0)
    end

    breaks = 20

    title_string = ht.title
    xlab_string = "Sample Size"
    ylab_string = string("Test power = 1 - β")
    legend_string4 = ""
    optimal_string3 = ""


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
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size d = ", d["d"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: Two-sample t test with n1 and n2
    elseif(ht.title == "T-test power calculation")
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
                return(NA)
            end
            return powerT2nTest(n1=na, n2=nb, d=d["d"], alpha = d["alpha"], alternative = d["alternative"])
        end

        power = [pwrt2test(ss) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size d = ", d["d"])
        legend_string3 = string("alpha = ", d["alpha"])
        legend_string4 = string("n1/n2 = ",round(n_rel,2))
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ",  d["n1"], " + ", d["n2"], " = ", n)

    # case: Difference of proportion (same sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ ss-> powerTwopTest(n=ss, h=d["h"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size h = ", d["h"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: difference of proportion (different sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
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
                return(NA)
            end
            return power2p2nTest(n1=na, n2=nb, h=d["h"], alpha = d["alpha"], alternative = d["alternative"])
        end

        power = [ss->pwr2p2ntest(ss) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size h = ", d["h"])
        legend_string3 = string("alpha = ", d["alpha"])
        legend_string4 = string("n1/n2 = ",round(n_rel,2))
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ",  d["n1"], " + ", d["n2"], " = ", n)

    # case: ANOVA
    elseif ht.title== "Balanced one-way analysis of variance power calculation"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerAnovaTest(n=ss, k=d["k"], f=d["f"], alpha = d["alpha"]) for ss in sample_sizes]

        # create labels
        legend_string1 = string("groups k = ", d["k"])
        legend_string2 = string("effect size f = ", d["f"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: Chi Squared
    elseif ht.title == "Chi squared power calculation"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerChisqTest(N=ss, w=d["w"], alpha = d["alpha"], df=d["df"]) for ss in sample_sizes]

        # create labels
        legend_string2 = string("effect size w = ", d["w"])
        legend_string1 = string("df = ", d["df"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("N = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: Normal distribution
    elseif ht.title == "Mean power calculation for normal distribution with known variance"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerNormTest(n=ss, d=d["d"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size d = ", d["d"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: proportion
    elseif ht.title == "Proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerPTest(n=ss, h=d["h"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size h = ", d["h"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
        optimal_string3 = d["note"]

    # case: correlation
    elseif ht.title == "Approximate correlation power calculation (arctangh transformation)"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerRTest(n=ss, r=d["r"], alpha = d["alpha"], alternative = d["alternative"]) for ss in sample_sizes]

        # create labels
        legend_string1 = string("tails = ", d["alternative"])
        legend_string2 = string("effect size r = ", d["r"])
        legend_string3 = string("alpha = ", d["alpha"])
        optimal_string1 = "optimal sample size"
        optimal_string2 = string("n = ", ceil(Int64,n))
    end

    # use DataFrame
    df = DataFrame(x=sample_sizes, y=power)
    df = df[completecases(df),:]
    len = size(df,1)

    # select the backend
    plotly()

    xlim_lower = df[len,:x] / 100

    # plot with title and x-axis and y-axis labels
    plot(df[:x],df[:y],
        marker=(2,.5,:circle,:blue),
        left_margin = 8mm,
        title = title_string,
        ylabel = ylab_string,
        xlabel = xlab_string,
        yticks=[0.0 0.2 0.4 0.6 0.8 1.0],
        ylims=(-0.03,1.0),
        xlims=(-(n_lower+xlim_lower),n_upper),
        label=false,
        legend=false,
        annotations=([((n_lower+xlim_lower),0.99,text(legend_string1,9,:blue,:left,:top)),
                ((n_lower+xlim_lower),0.94,text(legend_string2,9,:blue,:left,:top)),
                ((n_lower+xlim_lower),0.89,text(legend_string3,9,:blue,:left,:top)),
                ((n_lower+xlim_lower),0.84,text(legend_string4,9,:blue,:left,:top)),
            (sample_sizes[len-5],0.06,text(optimal_string1,9,:red,:left,:bottom)),
            (sample_sizes[len-5],0.01,text(string(optimal_string2,"  ", optimal_string3),9,:red,:left,:bottom))]))
    vline!([n],line=(1,:dot,.8,:orange))

end

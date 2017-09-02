using pwr, Plots, DataFrames

#f2 = pwr.F2Test(u=12,v=99,f2=.3,power=0.0)
import RecipesBase.plot

function plot(ht::pwr.htest)

    methods = (
        "One-sample t test power calculation",
        "Two-sample t test power calculation",
        "Paired t test power calculation",
        "T test power calculation",
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

    # case: One-sample, Two-sample or Paired t test
    if ht.title in ("One-sample t-test power calculation",
        "Two-sample t-test power calculation",
        "Paired t-test power calculation")

        n = d["n"]
        n_upper = ceil(Int64,max(n*1.5, n+30))
        n_increment = ceil(Int64,(n_upper - 10)/breaks)
        sample_sizes = collect(10:n_increment:n_upper)

        # obtain power at each sample size
        power = [ss->powerTTest(n=ss,d=d["d"],alpha = d["alpha"],sampletype=d["sampletype"],sided=d["alternative"]) for ss in sample_sizes]

        # create labels
        title_string = ht.title
        legend_string = string("tails =", d["alternative"], "\neffect size d =", d["d"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", ceil(Int64,n), "\n", d["note"])

    # case: Two-sample t test with n1 and n2
    elseif(ht.title == "T-test power calculation")
        n = d["n1"] + d["n2"]
        n_upper = ceil(Int64,max(n*1.5, n+30)) # upper at least 30 above n
        n_rel = d["n1"] / n # relative sample size; will be kept constant in claculations

        # generate data
        n_increment = ceil(Int64,(n_upper - 10)/breaks)
        sample_sizes = collect(10:n_increment:n_upper)

        function pwrt2test(ss)
            n1 = ceil(Int64,ss*n_rel)
            n2 = ss - n1
            if (n1 < 2 || n2 < 2)
                return(NA)
            end
            return powerT2nTest(n1=n1, n2=n2, d=d["d"], alpha = d["alpha"], sided = d["alternative"])
        end

        power = [ss->pwrt2test(ss) for ss in sample_sizes]

        # create labels
        title_string = ht.title
        legend_string = string("tails =", d["alternative"], "\neffect size d =", d["d"], "\nalpha =", d["alpha"], "\nn1/n2 = ", round(n_rel, 2))
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", d["n1"], " + ", d["n2"], " = ", n)

    # case: Difference of proportion (same sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ ss-> powerTwopTest(n=ss, h=d["h"], alpha = d["alpha"], sided = d["alternative"]) for ss in sample_sizes]

        # create labels
        title_string = "Difference of proportion power calculation\nfor binomial distribution (arcsine transformation)"
        legend_string = string("tails =", d["alternative"], "\neffect size h =", d["h"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", ceil(Int64,n), "\n", d["note"])

    # case: difference of proportion (different sample size)
    elseif ht.title == "Difference of proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["n1"] + d["n2"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_rel = d["n1"] / n # relative sample size; will be kept constant in claculations
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)

        function pwr2p2ntest(ss)
            n1 = ceil(Int64,ss*n_rel)
            n2 = ss - n1
            if (n1 < 2 || n2 < 2)
                return(NA)
            end
            return power2p2nTest(n1=n1, n2=n2, h=d["h"], alpha = d["alpha"], sided = d["alternative"])
        end

        power = [ss->pwr2p2ntest(ss) for ss in sample_sizes]

        # create labels
        title_string = "Difference of proportion power calculation\nfor binomial distribution (arcsine transformation)"
        legend_string = string("tails =", d["alternative"], "\neffect size h =", d["h"], "\nalpha =", d["alpha"], "\nn1/n2 = ", round(n_rel, 2))
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", d["n1"], " + ", d["n2"], " = ", n)

    # case: ANOVA
    elseif ht.title== "Balanced one-way analysis of variance power calculation"
        n = d["n"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerAnovaTest(n=ss, k=d["k"], f=d["f"], alpha = d["alpha"]) for ss in sample_sizes]

        # create labels
        title_string = "Balanced one-way analysis of variance \npower calculation"
        legend_string = string("groups k =", d["k"], "\neffect size f =", d["f"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string <- expression(paste("test power = 1 - ", beta))
        optimal_string <- paste("optimal sample size \nn = ", ceiling(n), "\n", x$note, sep = "")

    # case: Chi Squared
    elseif ht.title == "Chi squared power calculation"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerChisqTest(N=ss, w=d["w"], alpha = d["alpha"], df=d["df"]) for ss in sample_sizes]

        # create labels
        title_string = ht.title
        legend_string = string("effect size w =", d["w"], "\ndf =", d["df"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nN = ", ceil(Int64,n), "\n", d["note"])

    # case: Normal distribution
    elseif ht.title == "Mean power calculation for normal distribution with known variance"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerNormTest(n=ss, d=d["d"], alpha = d["alpha"], sided = d["alternative"]) for ss in sample_sizes]

        # create labels
        title_string = "Mean power calculation for normal distribution\nwith known variance"
        legend_string = string("tails =", d["alternative"], "\neffect size d =", d["d"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", ceil(Int64,n), "\n", d["note"])

    # case: proportion
    elseif ht.title == "Proportion power calculation for binomial distribution (arcsine transformation)"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerPTest(n=ss, h=d["h"], alpha = d["alpha"], sided = d["alternative"]) for ss in sample_sizes]

        # create labels
        title_string = "Proportion power calculation\nfor binomial distribution (arcsine transformation)"
        legend_string = string("tails =", d["alternative"], "\neffect size h =", d["h"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", ceil(Int64,n), "\n", d["note"])

    # case: correlation
    elseif ht.title == "approximate correlation power calculation (arctangh transformation)"
        n = d["N"]
        n_upper = ceil(Int64, max(n*1.5, n+30)) # upper at least 30 above n
        n_increment = ceil(Int64,(n_upper - 10)/breaks)

        # generate data
        sample_sizes = collect(10:n_increment:n_upper)
        power = [ss->powerRTest(n=ss, r=d["r"], alpha = d["alpha"], sided = d["alternative"]) for ss in sample_sizes]

        # create labels
        title_string = "approximate correlation power calculation\n(arctangh transformation)"
        legend_string = string("tails =", d["alternative"], "\nr =", d["r"], "\nalpha =", d["alpha"])
        xlab_string = "sample size"
        ylab_string = :(string("test power = 1 - ", beta))
        optimal_string = string("optimal sample size \nn = ", ceil(Int64,n))
    end

#  # pass arguments if required
#  if(length(dots <- list(...)) && !is.null(dots$xlab)){
#    xlab_string <- dots$xlab
#  }
#  if(length(dots <- list(...)) && !is.null(dots$ylab)){
#    ylab_string <- dots$ylab
#  }
#  if(length(dots <- list(...)) && !is.null(dots$main)){
#    title_string <- dots$main
# }
#
#  # position of text in plot
#  if(x$power < 0.5){
#    text_anchor <- 1
#    text_vjust <- 1
#  }else{
#    text_anchor <- 0
#    text_vjust <- 0
#  }
#  if(min(data$power, na.rm = TRUE) < 0.6){
#    legend_anchor <- 1
#    legend_vjust <- 1
#  }else{
#    legend_anchor <- 0
#    legend_vjust <- 0
#  }

    # use DataFrame
    df = DataFrame(x=sample_sizes, y=power)
    df = df[completecases(df),:]

    # select the backend
    gr()
    #plotly()

    # plot with title and x-axis and y-axis labels
    plot(df[:x],df[:y],title = title_string, xlabel = xlab_string,ylabel = ylab_string)

    # add options



end

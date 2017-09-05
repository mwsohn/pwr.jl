# test PowerPlot.jl
#

# T2nTest
using pwr, Roots
n1 = 220
n2 = 1832
n = n1 + n2
n_upper = ceil(Int64,max(n*1.5, n+30)) # upper at least 30 above n
n_rel = n1 / n # relative sample size; will be kept constant in claculations

# generate data
breaks = 20
n_increment = ceil(Int64,(n_upper - 10)/breaks)
n_lower = ceil(Int64,fzero(x->(x*n_rel) - 2.1,2,n2))
sample_sizes = collect(n_lower:n_increment:n_upper)

function pwrt2test(ss)
    na = ceil(Int64,ss*n_rel)
    nb = ss - na
    if (na < 2 || nb < 2)
        return(NA)
    end
    return powerT2nTest(n1=na, n2=nb, d=.2, alpha = .05, sided = "two")
end

power = [pwrt2test(ss) for ss in sample_sizes]

using Plots

legend_string1 = "tails = two sided"
legend_string2 = "effect size d = 0.2"
legend_string3 = "alpha = 0.05"
legend_string4 = string("n1/n2 = ",round(n_rel,2))
plot(sample_sizes,power,leg=false,
    title = "T-test power calculation",
    ylabel = "Test power = 1 - β",
    xlabel = "Sample Size",
    yticks=[0.0 0.2 0.4 0.6 0.8 1.0],
    ylim=(0.0,1.0),
    xlim=(10,n_upper),
    label=false,
    legend = false,
    annotations=([(50,0.99,text(legend_string1,9,:blue,:left,:top)),
            (50,0.94,text(legend_string2,9,:blue,:left,:top)),
            (50,0.89,text(legend_string3,9,:blue,:left,:top)),
            (50,0.84,text(legend_string4,9,:blue,:left,:top)),
        (sample_sizes[15],0.06,text(string("optimal sample size"),9,:red,:left,:bottom)),
        (sample_sizes[15],0.01,text(string("n = ",n1," + ",n2," = ",n),9,:red,:left,:bottom))]))
vline!([n],line=(1,:dot,.8,:red))




# plot with title and x-axis and y-axis labels
Plotly.plot(sample_sizes,power,
    title = "T-test power calculation",
    xlabel = "Test power = 1 - β",
    ylabel = "Sample Size",
    yticks=[0.0 0.2 0.4 0.6 0.8 1.0],
    ylim=(0.0,1.0),
    xlim=(10,n_upper),
    label=false,
    legend = false)


    ,
    annotations=([(20,0.99,text(legend_string,9,:blue,:left,:top)),
        (sample_sizes[15],0.05,text(optimal_string,9,:red,:left,:bottom))]))

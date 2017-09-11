function cohenES(;test::String = "",size::String = "")

    if test in ("p","t","r","anova","chisq","f2") == false
        error("`test` must be one of `p`, `t`, `r`, `anova`, `chisq`,or `f2`")
    end

    testd = Dict("p"=>1,"t"=>2,"r"=>3,"anova"=>4,"chisq"=>5,"f2"=>6)

    if in(size,("small","medium","large")) == false
        error("`size` must be one of `small`, `medium`, or `large`")
    end

    sized = Dict("small" => 1, "medium" => 2, "large" => 3)

    effsize = [[.2 .5 .8],
        [.2 .5 .8],
        [.1 .3 .5],
        [.1 .25 .4],
        [.1 .3 .5],
        [.02 .15 .35]]

    return htest(
        "Conventional effect size from Cohen (1982)",
        OrderedDict(
            "test" => test,
            "size" => size,
            "effectsize" => effsize[testd[test]][sized[size]]
            )
end

#cohenES("anova", "small")

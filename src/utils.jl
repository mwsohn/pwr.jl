using DataStructures

struct htest
    title::String
    d::OrderedDict
end

import Base.show

Base.show(io::IO, ::MIME"text/plain", h::htest) = show(io,h)

function show(io::IO, h::htest)
    print(io,"\n",h.title,"\n\n")

    for k in keys(h.d)
        if k == "note"
            continue
        end
        print(io,lpad(k,13)," = ",h.d[k],"\n")
    end

    if haskey(h.d,"note") && h.d["note"] != ""
        print(io,"\nNOTE: ", h.d["note"],"\n")
    end
end

function check_args(; kwargs...)
    for (key,val) in kwargs
        if key in (:h, :k, :n, :n1, :n2, :f, :f2, :r) && val == 0
            error("`",string(key),"` must be specified and greater than zero")
        elseif key in (:u, :v, :df)
            error("Degrees of freedom `",key,"` must be at least 1")
        elseif key == :w && val <= 0
            error("`w` must be positive")
        elseif key in (:alpha,:power) && (val <= 0.0 || val >= 1.0)
            error("`",key,"` must be a number in [0,1]")
        elseif key == :N && val < 1
            error("Number of observations `N` must be at least 1")
        end
    end
end

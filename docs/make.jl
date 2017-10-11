using Documenter, pwr

push!(LOAD_PATH,"../src/")

makedocs()

deplydocs(
    repo = "github.com/mwsohn/pwr.jl.git"
)

using Documenter, IERSConventions
using Pkg

const CI = get(ENV, "CI", "false") == "true"

if CI
    Pkg.add("Literate")
end

include("generate.jl")

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="IERSConventions.jl",
    modules=[IERSConventions],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "01 - Loading EOPs" => "Tutorials/gen/t01_load.md"
        ],
        "API" => [
            "Public API" => "API/api.md",
            "Low-level API" => "API/lapi.md"
        ]
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/IERSConventions.jl", branch="gh-pages"
)
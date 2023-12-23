using Documenter, IERS
using Pkg 

const CI = get(ENV, "CI", "false") == "true"


makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="IERS.jl",
    modules=[IERS],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "01 - Introduction" => "tutorials/t01_intro.md"
        ], 

        "Benchmarks" => "benchmarks.md",
        
        "API" => [
            "Public API" => "api/api.md",
            "Low-level API" => "api/lapi.md"
        ]
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
)
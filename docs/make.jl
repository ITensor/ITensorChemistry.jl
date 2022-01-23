using ITensorChemistry
using Documenter

DocMeta.setdocmeta!(ITensorChemistry, :DocTestSetup, :(using ITensorChemistry); recursive=true)

makedocs(;
    modules=[ITensorChemistry],
    authors="Matthew Fishman <mfishman@flatironinstitute.org> and contributors",
    repo="https://github.com/mtfishman/ITensorChemistry.jl/blob/{commit}{path}#{line}",
    sitename="ITensorChemistry.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mtfishman.github.io/ITensorChemistry.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mtfishman/ITensorChemistry.jl",
    devbranch="main",
)

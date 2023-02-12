using TransitionMatrices
using Documenter

DocMeta.setdocmeta!(TransitionMatrices, :DocTestSetup, :(using TransitionMatrices); recursive=true)

makedocs(;
    modules=[TransitionMatrices],
    authors="Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo="https://github.com/lucifer1004/TransitionMatrices.jl/blob/{commit}{path}#{line}",
    sitename="TransitionMatrices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucifer1004.github.io/TransitionMatrices.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucifer1004/TransitionMatrices.jl",
    devbranch="main",
)

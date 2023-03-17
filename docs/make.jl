using TransitionMatrices
using Documenter

DocMeta.setdocmeta!(TransitionMatrices, :DocTestSetup, :(using TransitionMatrices);
                    recursive = true)

makedocs(;
         modules = [TransitionMatrices],
         authors = "Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
         repo = "https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/blob/{commit}{path}#{line}",
         sitename = "TransitionMatrices.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://JuliaRemoteSensing.github.io/TransitionMatrices.jl",
                                  edit_link = "main",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "Usage" => "usage.md",
             "API" => "api.md",
         ])

deploydocs(;
           repo = "github.com/JuliaRemoteSensing/TransitionMatrices.jl",
           devbranch = "main")

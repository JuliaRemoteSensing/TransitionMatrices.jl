using TransitionMatrices
using Documenter

# ── Example Pluto notebooks ───────────────────────────────────────────────────
# The rendered pages (docs/src/examples/*.md, with plots embedded as base64 SVG)
# are committed and used as-is by default, so CI does not have to run the
# notebooks (which load Plots/GR and its native graphics stack — fragile on
# headless Linux runners). To refresh them, render locally with
#   BUILD_EXAMPLES=true julia --project=docs docs/make.jl
# and commit the regenerated docs/src/examples/*.md. Each notebook self-activates
# the examples/ environment in its own Pluto worker.
const NB_DIR = normpath(joinpath(@__DIR__, "..", "examples"))
const NB_OUT = joinpath(@__DIR__, "src", "examples")
const NOTEBOOKS = ["shapes_gallery.jl", "solver_landscape.jl",
    "angular_scattering.jl", "orientation_averaging.jl", "spectral_sensitivity.jl",
    "rain_radar.jl"]

function notebook_title(nb)
    in_markdown = false
    for line in eachline(joinpath(NB_DIR, nb))
        s = strip(line)
        if !in_markdown
            startswith(s, "md\"\"\"") || continue
            in_markdown = true
            s = strip(s[6:end])
        elseif startswith(s, "\"\"\"")
            in_markdown = false
            continue
        end
        startswith(s, "# ") && return String(strip(s[3:end]))
    end
    String(replace(splitext(nb)[1], "_" => " "))
end

const NB_PAGE_PATHS = ["examples/" * replace(nb, ".jl" => ".md") for nb in NOTEBOOKS]
const NB_PAGES = [
    notebook_title(nb) => path
    for (nb, path) in zip(NOTEBOOKS, NB_PAGE_PATHS)
]

function build_examples()
    mkpath(NB_OUT)
    bopts = BuildOptions(NB_DIR; output_format = documenter_output,
        use_distributed = false, previous_dir = NB_DIR)
    build_notebooks(bopts, NOTEBOOKS)
    for nb in NOTEBOOKS
        md = replace(nb, ".jl" => ".md")
        cp(joinpath(NB_DIR, md), joinpath(NB_OUT, md); force = true)
    end
end

if get(ENV, "BUILD_EXAMPLES", "false") == "true"
    # Only the render path needs PlutoStaticHTML/Pluto and a headless GR; the
    # default docs build (CI) uses the committed pages and loads none of this.
    ENV["GKSwstype"] = "100"
    using PlutoStaticHTML
    build_examples()
end

DocMeta.setdocmeta!(TransitionMatrices, :DocTestSetup, :(using TransitionMatrices);
    recursive = true)

makedocs(;
    modules = [TransitionMatrices],
    checkdocs = :exported,
    authors = "Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo = "https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/blob/{commit}{path}#{line}",
    sitename = "TransitionMatrices.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaRemoteSensing.github.io/TransitionMatrices.jl",
        edit_link = "main",
        size_threshold_ignore = NB_PAGE_PATHS,   # rendered notebooks embed large plots
        assets = String[]),
    pages = [
        "Home" => "index.md",
        "Theory & conventions" => "theory.md",
        "Usage" => "usage.md",
        "Examples" => vcat(["Overview" => "examples/index.md"], NB_PAGES),
        "Linearization" => "linearization.md",
        "Performance" => "performance.md",
        "API" => "api.md",
        "Methods & references" => "references.md"
    ])

if get(ENV, "DOCUMENTER_SKIP_DEPLOY", "false") != "true"
    deploydocs(;
        repo = "github.com/JuliaRemoteSensing/TransitionMatrices.jl",
        devbranch = "main")
end

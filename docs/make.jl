using Documenter, StableDQMC

makedocs(
    modules = [StableDQMC],
    doctest = false,
    sitename = "StableDQMC.jl",
    pages = [
        "Home" => "index.md",
    ],
    # assets = ["assets/custom.css", "assets/custom.js"]
)

deploydocs(
    repo = "github.com/crstnbr/StableDQMC.jl.git",
    push_preview = true,
    # target = "site",
)

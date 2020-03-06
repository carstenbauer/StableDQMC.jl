using Documenter, StableDQMC

makedocs(
    modules = [StableDQMC],
    doctest = false,
    sitename = "StableDQMC.jl",
    pages = [
        "Home" => "index.md",
        "Stable Inversion" => "inversions.md",
        "Matrix Decompositions" => "decompositions.md",
        "Exports" => ["exports/qr_udt.md", "exports/svd.md"],
    ],
    # assets = ["assets/custom.css", "assets/custom.js"]
)

deploydocs(
    repo = "github.com/crstnbr/StableDQMC.jl.git",
    push_preview = true,
    # target = "site",
)

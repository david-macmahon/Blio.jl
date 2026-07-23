using Documenter, Blio

# Load optional dependencies so that extension methods are loaded and their
# docstrings become available to Documenter.
using HDF5
using DataFrames
using DimensionalData

makedocs(
    sitename = "Blio.jl",
    authors = "David MacMahon",
    modules = [Blio],
    format = Documenter.HTML(
        canonical = "https://david-macmahon.github.io/Blio.jl/stable/",
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Filterbank" => "filterbank.md",
        "GuppiRaw" => "guppiraw.md",
        "Extensions" => "extensions.md",
    ],
    doctest = true,
    checkdocs = :exports,
    remotes = nothing,
)

deploydocs(
    repo = "github.com/david-macmahon/Blio.jl.git",
    devbranch = "main",
)

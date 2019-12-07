using Documenter
using BSplines

DocMeta.setdocmeta!(BSplines, :DocTestSetup, :(using BSplines))

makedocs(
         sitename = "BSplines.jl",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         modules = [BSplines],
         pages = [
                  "Home" => "index.md"
                  "Manual" => [
                               "The `BSplineBasis` type" => "basis.md"
                               "The `Spline` type" => "spline.md"
                               "Higher-level functions" => "functions.md"
                               "Plotting" => "plotting.md"
                              ]
                  "API documentation" => "api.md"
                 ]
        )

deploydocs(repo = "github.com/sostock/BSplines.jl.git")

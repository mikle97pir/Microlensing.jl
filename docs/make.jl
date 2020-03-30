using Documenter, Microlensing

makedocs(
    modules = [ Microlensing ],
    doctest = true,
    format = :html,
    sitename = "Microlensing.jl",
    authors = "Mikhail Pirogov",
    linkcheck = false,
    pages = [
        "Home" => "index.md"
    ],
    # Use clean URLs, unless built as a "local" build
    html_prettyurls = !("local" in ARGS),
)

deploydocs(
    repo    = "github.com/mikle97pir/mikle97pir.jl.git",
    target  = "build",
    deps    = nothing,
    make    = nothing
)

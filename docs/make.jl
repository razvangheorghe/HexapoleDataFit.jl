using Documenter, HexapoleDataFit

makedocs(;
    modules=[HexapoleDataFit],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/razvangheorghe/HexapoleDataFit.jl/blob/{commit}{path}#L{line}",
    sitename="HexapoleDataFit.jl",
    authors="Razvan Gheorghe, University of Oxford",
    assets=String[],
)

deploydocs(;
    repo="github.com/razvangheorghe/HexapoleDataFit.jl",
)

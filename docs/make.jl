using IgBLAST
using Documenter

DocMeta.setdocmeta!(IgBLAST, :DocTestSetup, :(using IgBLAST); recursive=true)

makedocs(;
    modules=[IgBLAST],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="IgBLAST.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/IgBLAST.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/IgBLAST.jl",
    devbranch="main",
)

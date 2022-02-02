using TemplateMatching
using Documenter

DocMeta.setdocmeta!(TemplateMatching, :DocTestSetup, :(using TemplateMatching); recursive=true)

makedocs(;
    modules=[TemplateMatching],
    authors="Stefano Campanella <15182642+stefanocampanella@users.noreply.github.com> and contributors",
    repo="https://github.com/stefanocampanella/TemplateMatching.jl/blob/{commit}{path}#{line}",
    sitename="TemplateMatching.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stefanocampanella.github.io/TemplateMatching.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stefanocampanella/TemplateMatching.jl",
    devbranch="master",
)

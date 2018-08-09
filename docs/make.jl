using Documenter, Fishmetica
 
makedocs(modules=[Fishmetica],
        doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/picestode/Fishmetica.git",
    julia  = "0.6.2",
    osname = "linux")
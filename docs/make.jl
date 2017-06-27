using Documenter, LowFrequencyEconometrics

makedocs(modules=LowFrequencyEconometrics,
        doctest=true)

deploydocs(deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "https://github.com/p-chaim/LowFrequencyEconometrics.jl.git",
    julia  = "0.5.2",
    osname = "windows")

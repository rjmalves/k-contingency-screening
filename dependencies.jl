using Pkg

dependencies = [
    "Distributed",
    "SharedArrays",
    "LightGraphs",
    "GraphIO",
    "Printf",
    "LinearAlgebra",
    "Combinatorics",
    "IterTools",
    "ArgParse",
    "StatsBase",
    "GraphPlot"
]

Pkg.add(dependencies)

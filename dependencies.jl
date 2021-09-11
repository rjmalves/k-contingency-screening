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
    "GraphPlot",
    "BenchmarkTools"
]

Pkg.add(dependencies)

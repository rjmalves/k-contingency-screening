using LightGraphs, GraphIO
include("cfb.jl")
using Main.CFB

g = Graph(loadgraph("t39.txt", "Grafo", EdgeListFormat()))

using GraphPlot
gplot(g, nodelabel=1:nv(g))

cfb_guloso(g, 3)

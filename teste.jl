using LightGraphs, GraphIO
include("cfb.jl")
using Main.CFB

g = Graph(loadgraph("t11.txt", "Grafo", EdgeListFormat()))

r, b, m = iteracao_cfb(g)


using GraphPlot
gplot(g, nodelabel=1:nv(g))

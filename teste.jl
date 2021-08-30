using LightGraphs, GraphIO
include("cfb.jl")
using Main.CFB

g = Graph(loadgraph("t11.txt", "Grafo", EdgeListFormat()))

m = iteracao_cfb(g)


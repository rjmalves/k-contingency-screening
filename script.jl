include("cfb.jl")
using Main.CFB, BenchmarkTools

println("------ ieee300 -------")
g = read_edgelist("ieee300.txt")
@time cfb_exaustivo(g, 1, "ieee300")
@time cfb_exaustivo(g, 2, "ieee300")
@time cfb_exaustivo(g, 3, "ieee300")
@time cfb_exaustivo(g, 4, "ieee300")

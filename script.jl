include("cfb.jl")
using Main.CFB, BenchmarkTools

println("------ ieee118 -------")
g = read_edgelist("ieee118.txt")
@time cfb_exaustivo(g, 1, "ieee118")
# @time cfb_exaustivo(g, 2, "ieee118")
# @time cfb_exaustivo(g, 3, "ieee118")
# @time cfb_exaustivo(g, 4, "ieee118")

# println("------ t300 -------")
# g = read_edgelist("t300.txt")
# @time cfb_exaustivo(g, 1, "t300")
# @time cfb_exaustivo(g, 2, "t300")
# @time cfb_exaustivo(g, 3, "t300")
# @time cfb_exaustivo(g, 4, "t300")

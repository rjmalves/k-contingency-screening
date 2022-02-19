include("cfb.jl")
using Main.CFB, BenchmarkTools


println("------ t11 -------")
g = read_edgelist("t11.txt")
@time cfb_exaustivo(g, 1, "t11")
@time cfb_exaustivo(g, 2, "t11")
@time cfb_exaustivo(g, 3, "t11")
@time cfb_exaustivo(g, 4, "t11")

println("------ t39 -------")
g = read_edgelist("t39.txt")
@time cfb_exaustivo(g, 1, "t39")
@time cfb_exaustivo(g, 2, "t39")
@time cfb_exaustivo(g, 3, "t39")
@time cfb_exaustivo(g, 4, "t39")

println("------- t57 --------")
g = read_edgelist("t57.txt")
@time cfb_exaustivo(g, 1, "t57")
@time cfb_exaustivo(g, 2, "t57")
@time cfb_exaustivo(g, 3, "t57")
@time cfb_exaustivo(g, 4, "t57")

println("------ t118 -------")
g = read_edgelist("t118.txt")
@time cfb_exaustivo(g, 1, "t118")
@time cfb_exaustivo(g, 2, "t118")
@time cfb_exaustivo(g, 3, "t118")
@time cfb_exaustivo(g, 4, "t118")

println("------ t300 -------")
g = read_edgelist("t300.txt")
@time cfb_exaustivo(g, 1, "t300")
@time cfb_exaustivo(g, 2, "t300")
@time cfb_exaustivo(g, 3, "t300")
@time cfb_exaustivo(g, 4, "t300")

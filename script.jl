include("cfb/cfb.jl")
using Main.CFB, BenchmarkTools

# Escolha do arquivo a ser utilizado
sistema = "itaipu11"
arquivo = "$sistema.txt"


# Chamadas ao cálculo da CFB
println("------ $sistema -------")
g = read_edgelist(arquivo)
@time cfb_exaustivo(g, 1, sistema)
@time cfb_exaustivo(g, 2, sistema)
@time cfb_exaustivo(g, 3, sistema)
@time cfb_exaustivo(g, 4, sistema)

using Graphs, ArgParse, BenchmarkTools


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "grafo"
            help = "Arquivo com o grafo em lista de arestas"
            arg_type = String
            required = true
        "k"
            help = "Número de arestas removidas simultâneamente"
            arg_type = Int
            required = true
    end

    return parse_args(s)
end

# # Echo das entradas
# args = parse_commandline()
# for (arg, val) in args
#     println("$arg  =>  $val")
# end

# ARQ = args["grafo"]
# K = args["k"]
# NOME = string(split(ARQ, ".")[1])

ARQ = "t39.txt"
K = 1
NOME = string(split(ARQ, ".")[1])


include("cfb.jl")
using Main.CFB

using GraphPlot

g = read_edgelist(ARQ)
cfb_exaustivo(g, K, NOME)

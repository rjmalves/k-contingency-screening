using LightGraphs, GraphIO, ArgParse, BenchmarkTools


# function parse_commandline()
#     s = ArgParseSettings()
#     @add_arg_table s begin
#         "grafo"
#             help = "Arquivo com o grafo em lista de arestas"
#             arg_type = String
#             required = true
#         "k"
#             help = "Número de arestas removidas simultâneamente"
#             arg_type = Int
#             required = true
#     end

#     return parse_args(s)
# end

# # Echo das entradas
# args = parse_commandline()
# for (arg, val) in args
#     println("$arg  =>  $val")
# end

# ARQ = args["grafo"]
# K = args["k"]

include("cfb.jl")
using Main.CFB

using GraphPlot

for ARQ in ["ieee118.txt"]
    for K in [2, 3, 4]
        println(ARQ, " ", K)
        NOME = string(split(ARQ, ".")[1])
        g = Graph(loadgraph(ARQ, NOME, EdgeListFormat()))
        @btime cfb_exaustivo($g, $K, $NOME)
    end
end
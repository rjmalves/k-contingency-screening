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

# Echo das entradas
args = parse_commandline()
for (arg, val) in args
    println("$arg  =>  $val")
end

ARQ = args["grafo"]
K = args["k"]

include("cfb.jl")
using Main.CFB

NOME = string(split(ARQ, ".")[1])
g = read_edgelist(ARQ)
cfb_guloso(g, K, NOME)

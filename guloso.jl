using Distributed, LightGraphs, GraphIO, Printf
using SharedArrays, IterTools, Combinatorics, ArgParse


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chunk_lim"
           help = "Memória disponível para execução (MB)"
           arg_type = Int
           default = 1000
        "grafo"
            help = "Arquivo com o grafo em lista de arestas"
            arg_type = String
            required = true
        "nprocs"
            help = "Número de processos paralelos"
            arg_type = Int
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
NPROCS = args["nprocs"]
K = args["k"]
chunk = args["chunk_lim"]

NOME = string(split(ARQ, ".")[1])
addprocs(NPROCS)

@everywhere include("cfb.jl")
@everywhere using Main.CFB

g = Graph(loadgraph(ARQ, NOME, EdgeListFormat()))
cfb_guloso(g, K, chunk, NOME)

using Distributed, LightGraphs, GraphIO, Printf
using SharedArrays, IterTools, Combinatorics, ArgParse


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--iter_num"
           help = "Número de iterações para execução do DE"
           arg_type = Int
           default = 10
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
        "pop"
            help = "Tamanho da população do DE"
            arg_type = Int
            required = true
        "crossover"
            help = "Taxa de crossover do DE"
            arg_type = Float64
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
iter_num = args["iter_num"]
pop = args["pop"]
crossover = args["crossover"]

NOME = string(split(ARQ, ".")[1])
addprocs(NPROCS)

@everywhere include("cfb.jl")
@everywhere using Main.CFB

g = Graph(loadgraph(ARQ, NOME, EdgeListFormat()))
cfb_de(g, K, pop, crossover, iter_num)

using Graphs, ArgParse, BenchmarkTools


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--iter_num"
           help = "Número de iterações para execução do DE"
           arg_type = Int
           default = 10
        "--beta_min"
           help = "Valor mínimo do parâmetro de mutação"
           arg_type = Float64
           default = 0.2
        "--beta_max"
           help = "Valor máximo do parâmetro de mutação"
           arg_type = Float64
           default = 0.8
        "grafo"
            help = "Arquivo com o grafo em lista de arestas"
            arg_type = String
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
K = args["k"]
iter_num = args["iter_num"]
beta_min = args["beta_min"]
beta_max = args["beta_max"]
pop = args["pop"]
crossover = args["crossover"]


include("cfb.jl")
using Main.CFB

NOME = string(split(ARQ, ".")[1])
g = read_edgelist(ARQ)
cfb_de(g, K, pop, crossover, beta_min, beta_max, iter_num, NOME)

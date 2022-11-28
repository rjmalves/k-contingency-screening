module CFB

using Graphs: Graph, SimpleGraph, nv, ne, adjacency_matrix
using Graphs: incidence_matrix, loadgraph, savegraph
using Graphs: rem_edge!, add_edge!, add_vertices!
using Graphs: is_connected, edges, src, dst
using Graphs: betweenness_centrality, closeness_centrality
using LinearAlgebra: transpose, diagm, inv
using IterTools: subsets
using StatsBase: sample, sortperm
using Combinatorics: multinomial
using CSV, DataFrames, Tables
using Random: seed!

export read_edgelist
export cfb_exaustivo
export cfb_de
export cfb_guloso
export current_flow_betweenness

# Parâmetros: grafo, k
seed!(2021)


function read_edgelist(s::String)
    f = open(s, "r")

    g = SimpleGraph(Int64)

    nodes = Set([])
    edges = Set([])

    while !eof(f)
        line = readline(f)
        if contains(line, " ")
            char = " "
        elseif contains(line, "\t")
            char = "\t"
        end
        src = -1
        dst = -1
        for node = split(line, char)
            if length(node) == 0
                continue
            end
            n = strip(node)
            if src == -1
                src = parse(Int64, n)
            else
                dst = parse(Int64, n)
            end
        end
        push!(nodes, src)
        push!(nodes, dst)
        push!(edges, [src, dst])
    end
    close(f)

    add_vertices!(g, length(nodes))
    for e = edges
        add_edge!(g, e[1], e[2])
    end

   return g
end


function current_flow_betweenness(g::Graph)
    # Calcula a current flow betweenness segundo o algoritmo de Brandes
    n = nv(g)
    m = ne(g)
    betweenness = zeros(Float64, n)
    # A matriz de adjacência
    A = adjacency_matrix(g)
    # A matriz de incidência no formato a ser utilizado
    b = -transpose(incidence_matrix(g, oriented=true))
    # Constante de normalização
    n_b = (n - 1) * (n - 2)
    # Tomando a laplaciana do grafo, invertendo sua submatriz
    # e preenchendo com zeros
    L = diagm(0 => reduce(+, A, dims=[1])[1:n]) - A
    L_tilde = L[2:n, 2:n]
    C = vcat(zeros(1, n), hcat(zeros(n-1, 1), inv(L_tilde)))
    # A matriz de fluxo de Brandes
    F = b*C
    # Iterando na matriz de fluxo (linha a linha)
    for j = 1:m
        # Ordena as entradas da linha em ordem
        # não-decrescente e pega os índices dos índices
        pos = sortperm(sortperm(-F[j, :]))
        # Encontra o índice da fonte e do dreno,
        # segundo a incidência
        e = findall(b[j, :] .!= 0)
        s = e[1]
        d = e[2]
        # Atualiza a betweenness para cada vértice
        for i = 1:n
            betweenness[s] += (i - pos[i]) * F[j, i]
            betweenness[d] += (n + 1 - i - pos[i]) * F[j, i]
        end
    end
    # Normaliza a betweenness
    norm_vec = 1:n
    betweenness -= norm_vec
    betweenness .+= 1
    betweenness .*= 2 / n_b
    return betweenness
end

# function current_flow_betweenness(g::Graph)
#     return betweenness_centrality(g)
# end

function edge_indices_k_tuples(m::Integer, k::Integer)
    subs = subsets(1:m, k)
    edge_index_matrix = zeros(Integer, length(subs), k)
    for (i, s) = enumerate(subs)
        edge_index_matrix[i, :] = s
    end
    return edge_index_matrix
end

function edges_k_tuples(g::Graph, indices::Matrix{Integer})
    edgs = [e for e = edges(g)]
    n_subsets, size_subsets = size(indices)
    edge_matrix = zeros(Integer, n_subsets, 2 * size_subsets)
    for i = 1:n_subsets
        for (j, e) = enumerate(indices[i, :])
            edg = edgs[e]
            edge_matrix[i, 2 * j - 1] = src(edg)
            edge_matrix[i, 2 * j] = dst(edg)
        end
    end
    return edge_matrix
end

function edit_graph(g::Graph, edgs::Vector{Integer})
    g_bkp = Graph(g)
    n_edges = length(edgs)
    n_edges = Int64(n_edges / 2)
    for i = 1:n_edges
        rem_edge!(g_bkp, edgs[2 * i - 1], edgs[2 * i])
    end
    return g_bkp
end

function check_valid_removals(g::Graph, tuples::Matrix{Integer})
    n_subsets, size_subsets = size(tuples)
    edges_by_subset = Int64(size_subsets / 2)
    valids = [true for i = 1:n_subsets]
    for i = 1:n_subsets
        g_edit = edit_graph(g, tuples[i, :])
        valids[i] = is_connected(g_edit)
    end
    return valids
end

function filter_valid_removals(tuples::Matrix{Integer}, valids::Vector{Bool})
    n_subsets, size_subsets = size(tuples)
    valid_indices = [i for i = 1:n_subsets][valids]
    num_valids = Int64(sum(valids))
    edge_tuples_matrix = zeros(Integer, num_valids, size_subsets)
    for (i, ti) = enumerate(valid_indices)
        edge_tuples_matrix[i, :] = tuples[ti, :]
    end
    return edge_tuples_matrix
end

function filter_invalid_removals(tuples::Matrix{Integer}, valids::Vector{Bool})
    n_subsets, size_subsets = size(tuples)
    invalids = .!valids
    invalid_indices = [i for i = 1:n_subsets][invalids]
    num_invalids = Int64(sum(invalids))
    edge_tuples_matrix = zeros(Integer, num_invalids, size_subsets)
    for (i, ti) = enumerate(invalid_indices)
        edge_tuples_matrix[i, :] = tuples[ti, :]
    end
    return edge_tuples_matrix
end

function betweenness_from_removals(g::Graph, tuples::Matrix{Integer})
    n = nv(g)
    n_subsets, _ = size(tuples)
    betweenness = zeros(Float64, n_subsets, n)
    for i = 1:n_subsets
        g_edit = edit_graph(g, tuples[i, :])
        betweenness[i, :] = current_flow_betweenness(g_edit)
    end
    return betweenness
end

function betweenness_deltas(betweenness::Matrix{Float64},
                            reference::Vector{Float64})
    n_subsets, n = size(betweenness)
    deltas = zeros(Float64, n_subsets, n)
    for i = 1:n_subsets
        deltas[i, :] = abs.(betweenness[i, :] - reference)
    end
    return deltas
end

function deltas_global_by_removal(betweenness::Matrix{Float64})
    return abs.(sum(betweenness, dims=2))
end

function normalized_deltas_global_by_vertex(deltas::Matrix{Float64})
    global_deltas = sum(deltas, dims=1)
    removal_bets = abs.(global_deltas)
    # minimo = minimum(removal_bets[removal_bets .> 1e-10])
    # removal_bets = (removal_bets .- minimo) ./ (maximum(removal_bets) - minimo)
    normalized = vec(removal_bets)
    return DataFrame(V=1:length(normalized), DELTA=normalized)
end

function normalized_deltas_global_by_edge(deltas::Matrix{Float64},
                                          edges::Matrix{Integer},
                                          tuples::Matrix{Integer})
    n_subsets, n = size(deltas)
    _, k = size(tuples)
    m, _ = size(edges)
    k = Int64(k / 2)
    deltas_by_tuple = sum(abs.(deltas), dims=2)
    deltas_by_tuple = abs.(deltas_by_tuple)
    deltas_by_edge = zeros(Float64, m)
    # Makes a map between the edge vertex pair and the index
    d = Dict([0,0]=>0)
    for i = 1:m
        d[edges[i, :]] = i
    end
    for i = 1:n_subsets
        for j = 1:k
            edg = [tuples[i, 2 * j - 1], tuples[i, 2 * j] ]
            idx = d[edg]
            deltas_by_edge[idx] += deltas_by_tuple[i]
        end
    end
    # deltas_by_edge = deltas_by_edge .- minimum(deltas_by_edge)
    # deltas_by_edge = deltas_by_edge ./ (maximum(deltas_by_edge) - minimum(deltas_by_edge))
    df = DataFrame(SRC=edges[:, 1], DST=edges[:, 2], DELTA=vec(deltas_by_edge))
    return df
end

function normalized_deltas_local(deltas::Matrix{Float64})
    n_subsets, n = size(deltas)
    normalized = abs.(sum(deltas, dims=2))
    # normalized = normalized .- minimum(normalized)
    # normalized = normalized ./ maximum(normalized)
    return vec(normalized)
end

function export_exaustive_results(g::Graph,
                                  valid_tuples::Matrix{Integer},
                                  globals::DataFrame,
                                  edge_globals::DataFrame,
                                  locals::Vector{Float64},
                                  disconnects::Matrix{Integer},
                                  k::Integer,
                                  filename::String,
                                  metodo::String)
    dir = string(metodo, "_", filename, "_", string(k))
    dir_bkp = pwd()
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    CSV.write("valid_tuples.csv", Tables.table(valid_tuples), writeheader=false)
    CSV.write("disconnects.csv", Tables.table(disconnects), writeheader=false)
    CSV.write("local_deltas.csv", Tables.table(locals), writeheader=false)
    CSV.write("vertex_global_deltas.csv", globals, writeheader=false)
    CSV.write("edge_global_deltas.csv", edge_globals, writeheader=false)
    cd(dir_bkp)
end

function export_de_results(g::Graph,
                           valid_tuples::Matrix{Integer},
                           globals::Vector{Float64},
                           edge_globals::Vector{Float64},
                           locals::Vector{Float64},
                           disconnects::Matrix{Integer},
                           k::Integer,
                           filename::String,
                           metodo::String)
    dir = string(metodo, "_", filename, "_", string(k))
    dir_bkp = pwd()
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    savegraph("edgelist.csv", g, EdgeListFormat())
    CSV.write("valid_tuples.csv", Tables.table(valid_tuples), writeheader=false)
    CSV.write("disconnects.csv", Tables.table(disconnects), writeheader=false)
    CSV.write("local_deltas.csv", Tables.table(locals), writeheader=false)
    CSV.write("vertex_global_deltas.csv", Tables.table(globals), writeheader=false)
    CSV.write("edge_global_deltas.csv", Tables.table(edge_globals), writeheader=false)
    cd(dir_bkp)
end

function export_greedy_results(g::Graph,
                               edges::Matrix{Integer},
                               deltas::Vector{Float64},
                               filename::String)
    dir = string("guloso", "_", filename)
    dir_bkp = pwd()
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    savegraph("edgelist.csv", g, EdgeListFormat())
    CSV.write("edges.csv", Tables.table(edges), writeheader=false)
    CSV.write("global_deltas.csv", Tables.table(deltas), writeheader=false)
    cd(dir_bkp)
end

function cfb_exaustivo(g::Graph, k::Integer, arquivo_saida="result"::String)
    # 1 - Calcula a betweenness de referência
    ref_cfb = current_flow_betweenness(g)
    # 2 - Obtém as tuplas de k arestas
    n = nv(g)
    m = ne(g)
    edgs = edges_k_tuples(g, edge_indices_k_tuples(m, 1))
    s = edge_indices_k_tuples(m, k)
    e = edges_k_tuples(g, s)
    # 3 - Verifica quais remoções desconectam o grafo - salva como
    #     cenários mais críticos e não avaliados 
    v = check_valid_removals(g, e)
    ve = filter_valid_removals(e, v)
    ive = filter_invalid_removals(e, v)
    # 4 - Para cada remoção que não desconecta o grafo, calcula a betweenness
    bets = betweenness_from_removals(g, ve)
    # 5 - Calcula os deltas da betweenness
    deltas = betweenness_deltas(bets, ref_cfb)
    # 7 - Para cada barra, realiza a soma dos deltas e normaliza
    #     para o intervalo 0 - 1 (impacto global)
    global_norms = normalized_deltas_global_by_vertex(deltas)
    global_edge_norms = normalized_deltas_global_by_edge(deltas, edgs, ve)
    local_norms = normalized_deltas_local(deltas)
    # 8 - Escreve em arquivo de texto a matriz de deltas
    export_exaustive_results(g, ve, global_norms, global_edge_norms,
                             local_norms, ive, k, arquivo_saida, "exaustivo")
end

function sample_initial_pop(m::Integer,
                            n_pop::Integer,
                            k::Integer)
    indices_pop = zeros(Integer, n_pop, k)
    indices = 1:m
    for i = 1:n_pop
        indices_pop[i, :] = sample(indices, k, replace=false)    
    end
    return indices_pop
end

function de_cost_function(deltas::Matrix{Float64})
    n_subsets, n = size(deltas)
    costs = zeros(Float64, n_subsets, n)
    for i = 1:n_subsets
        removal_bets = abs.(deltas[i, :])
        costs[i, :] = removal_bets
    end
    return vec(sum(costs, dims=2))
end

function de_mutation(m::Integer,
                     k::Integer,
                     pop_size::Integer,
                     indices_pop::Matrix{Integer},
                     beta_min=0.2::Float64,
                     beta_max=0.8::Float64)
    coefs = rand(-m:m, pop_size, k)
    beta = rand(pop_size, k) .* (beta_max - beta_min) .+ beta_min
    indices_pop_mut = round.(Integer,
                             indices_pop .+ (beta .* coefs))
    indices_pop_mut = max.(ones(Integer, pop_size, k),
                           indices_pop_mut)
    indices_pop_mut = min.(m * ones(Integer, pop_size, k),
                           indices_pop_mut)
    return indices_pop_mut
end

function de_crossover(pop_size::Integer,
                      pop_indices::Matrix{Integer},
                      pop_indices_mut::Matrix{Int64},
                      crossover_rate=0.5::Float64)
    flags = rand(pop_size, 1) .< crossover_rate
    for c = 1:pop_size
        if flags[c]
            pop_indices[c, :] = pop_indices_mut[c, :]
        end
    end
    return pop_indices
end

function de_iter!(g::Graph,
                  pop_indices::Matrix{Integer},
                  ref_bets::Vector{Float64},
                  crossover_rate::Float64,
                  beta_min::Float64,
                  beta_max::Float64,
                  ultima=false::Bool
                  )
    m = ne(g)
    # 1 - Extrai as tuplas de k arestas
    edgs = edges_k_tuples(g, edge_indices_k_tuples(m, 1))
    e = edges_k_tuples(g, pop_indices)
    n_pop, k = size(pop_indices)
    # 2 - Filtra as tuplas que não desconectam o grafo
    v = check_valid_removals(g, e)
    while sum(v) == 0
        pop_indices = sample_initial_pop(m, n_pop, k)
        e = edges_k_tuples(g, pop_indices)
        v = check_valid_removals(g, e)
    end
    ve = filter_valid_removals(e, v)
    ive = filter_invalid_removals(e, v)
    # 3 - Calcula as betweenness das remoções válidas
    bets = betweenness_from_removals(g, ve)
    # 4 - Calcula os deltas da betweenness
    deltas = betweenness_deltas(bets, ref_bets)
    # 5 - Calcula a função de custo do DE
    costs = de_cost_function(deltas)
    # 6 - Ordena as tuplas por custo
    p = sortperm(-costs)
    costs = costs[p]
    pop_indices = pop_indices[p, :]
    pop_size, k = size(pop_indices)
    if !ultima
        # 7 - Elitismo
        best_edg = pop_indices[1, :]
        # 8 - Mutação
        pop_indices_mut = de_mutation(m,
                                      k,
                                      pop_size,
                                      pop_indices,
                                      beta_min,
                                      beta_max)
        # 9 - Crossover
        pop_indices_cross = de_crossover(pop_size,
                                         pop_indices,
                                         pop_indices_mut,
                                         crossover_rate)
        # 10 - Atualiza a População
        pop_indices = pop_indices_cross
        pop_indices[1, :] = best_edg
        return nothing, nothing, nothing, nothing
    else
        # Calcula os deltas da população final
        global_norms = normalized_deltas_global_by_vertex(deltas)
        global_edge_norms = normalized_deltas_global_by_edge(deltas, edgs, ve)
        local_norms = normalized_deltas_local(deltas)
        return ve, ive, global_norms, global_edge_norms, local_norms
    end
end

function cfb_de(g::Graph,
                k::Integer,
                pop_size::Integer,
                crossover_rate::Float64,
                beta_min::Float64,
                beta_max::Float64,
                iter_num=10::Integer,
                arquivo_saida="result"::String)
    # 1 - Calcula a betweenness de referência
    ref_cfb = current_flow_betweenness(g)
    # 2 - Obtém as tuplas de índices das k arestas
    m = ne(g)
    edge_index_pop = sample_initial_pop(m, pop_size, k)
    for i = 1:iter_num-1
        de_iter!(g, edge_index_pop, ref_cfb,
                 crossover_rate, beta_min, beta_max)
    end
    # A última retorna os impactos normalizados da população
    ve, ive, global_norms, edge_norms, local_norms = de_iter!(g,
                                                              edge_index_pop,
                                                              ref_cfb,
                                                              crossover_rate,
                                                              beta_min,
                                                              beta_max,
                                                              true)
    # Escreve em arquivo de texto a matriz de deltas
    export_de_results(g, ve, global_norms, edge_norms,
                      local_norms, ive, k, arquivo_saida, "de")
end

function iteracao_cfb_guloso(g::Graph)
    # 1 - Calcula a betweenness de referência
    ref_cfb = current_flow_betweenness(g)
    # 2 - Obtém as tuplas de k = 1 arestas
    k = 1
    n = nv(g)
    m = ne(g)
    s = edge_indices_k_tuples(m, k)
    e = edges_k_tuples(g, s)
    # 3 - Verifica quais remoções desconectam o grafo - salva como
    #     cenários mais críticos e não avaliados 
    v = check_valid_removals(g, e)
    ve = filter_valid_removals(e, v)
    ive = filter_invalid_removals(e, v)
    # 4 - Para cada remoção que não desconecta o grafo, calcula a betweenness
    bets = betweenness_from_removals(g, ve)
    # 5 - Calcula os deltas da betweenness
    deltas = betweenness_deltas(bets, ref_cfb)
    # 7 - Para cada barra, realiza a soma dos deltas e normaliza
    #     para o intervalo 0 - 1 (impacto global)
    global_deltas = deltas_global_by_removal(deltas)
    max_delta, idx = findmax(vec(global_deltas))
    return ve[idx, :], max_delta
end

function cfb_guloso(g::Graph,
                    k::Integer,
                    arquivo_saida="result"::String)
    g_bkp = Graph(g)
    edges = zeros(Integer, k, 2)
    deltas = zeros(Float64, k)
    for i = 1:k
        edg, delta = iteracao_cfb_guloso(g_bkp)
        nodelabel = 1:nv(g_bkp)
        edges[i, :] = edg
        deltas[i] = delta
        rem_edge!(g_bkp, edg[1], edg[2])
    end
    export_greedy_results(g, edges, deltas, arquivo_saida)
end

end
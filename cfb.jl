module CFB

using Distributed, SharedArrays, LightGraphs
using GraphIO, Printf, LinearAlgebra, IterTools
using StatsBase, Combinatorics

export cfb_exaus
export cfb_de
export iteracao_cfb


function edit_adjacencies!(adjacencies::SharedArray,
                           A_temp,
                           edg_tuples,
                           N::Int64,
                           calc_range)
    println("Editando adjacências")
    # Função interna da paralelização para edição da Matriz
    # de adjacência por remoção de arestas
    i_idx = calc_range[1]
    f_idx = calc_range[2]
    # Para o intervalo de arestas recebido pelo processo
    for idx = i_idx:f_idx
        # Inicia com a matriz de adjacência original
        adjacencies[:, :, idx] = A_temp
        # Para cada aresta a ser removida, zera as entradas da matriz
        for eidx = 1:N
            i, j = edg_tuples[idx][eidx][1], edg_tuples[idx][eidx][2]      
            adjacencies[i, j, idx] = 0
            adjacencies[j, i, idx] = 0
        end
    end
end

function verify_connectivity!(adjacencies::SharedArray,
                              tuples_to_remove::SharedArray,
                              calc_range)
    println("Verificando conectividade")
    # Função interna da paralelização para verificar
    # se os grafos editados são conexos
    i_idx = calc_range[1]
    f_idx = calc_range[2]
    # Para o intervalo de matrizes recebido pelo processo
    for idx = i_idx:f_idx
        # Se for conexo, sinaliza no vetor
        if is_connected(Graph(adjacencies[:, :, idx]))
            tuples_to_remove[idx] = true
        else
            # Senão, mantém falso
            tuples_to_remove[idx] = false
        end
    end
end

function current_flow_betweenness!(adjacencies::SharedArray,
                                   betweenness::SharedArray,
                                   calc_range)
    println("Calculando CFB")
    # Calcula a current flow betweenness segundo o algoritmo de Brandes
    i_idx = calc_range[1]
    f_idx = calc_range[2]
    # Para o intervalo de matrizes recebido pelo processo
    for idx = i_idx:f_idx
        # A matriz de adjacência
        A = adjacencies[:, :, idx]
        # A matriz de incidência no formato a ser utilizado
        b = -transpose(incidence_matrix(Graph(A), oriented=true))
        # Constantes importantes (número de nós, arestas
        # e constante de normalização)
        n = size(A)[1]
        m = size(b)[1]
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
                betweenness[s, idx] += (i - pos[i]) * F[j, i]
                betweenness[d, idx] += (n + 1 - i - pos[i]) * F[j, i]
            end
        end
        # Normaliza a betweenness
        for i in 1:n
            betweenness[i, idx] = (betweenness[i, idx] - i + 1) * 2 / n_b
        end
    end
end

"""
define_critical_nodes(noncritical_nodes,
                      medium_nodes,
                      critical_nodes,
                      betweenness,
                      reference_betweenness,
                      calc_range)

Realiza a comparação das diferenças das CFB em % para cada
configuração por remoção de arestas e categoriza cada nó como
crítico, médio ou não crítico.

### Argumentos

- `betweenness` que possui dimensão NxVE, onde N é o número
de nós do grafo e VE é o número de tuplas de arestas cujas remoções
são válidas, isto é, o grafo após a remoção permanece conexo. Este
array deve conter, na posição [i, j] o valor da betweenness do nó
i na configuração por remoção da j-ésima tupla de arestas válida.

- `reference_betweenness` deve possuir dimensão Nx1, onde N é o
número de nós do grafo. A posição [i] deve conter a CFB do nó i
no grafo original.

E retorna as informações através de escritas nas arrays:
`noncritical_nodes`, `medium_nodes` e `critical_nodes`, que são
de dimensão NxVE, onde N é o número de nós do grafo e VE é o número
de tuplas válidas para remoção. Se o nó i for sinalizado, respectivamente,
como não crítico (Delta CFB < 20%), médio (20% <= DeltaCFB <= 70%) ou crítico
(Delta CFB > 70%), mediante a remoção da tupla j, então a respectiva
matriz terá valor 1 na posição [i, j].

"""
function define_critical_nodes!(noncritical_nodes::SharedArray,
                                medium_nodes::SharedArray,
                                critical_nodes::SharedArray,
                                betweenness::SharedArray,
                                reference_betweenness::SharedArray,
                                calc_range)
    println("Definindo nós críticos")
    # Função interna da paralelização para obtenção
    # dos nós críticos em uma remoção
    i_idx = calc_range[1]
    f_idx = calc_range[2]
    # Descobre o número de nós
    n = length(reference_betweenness)
    # Para o intervalo de remoções recebido pelo processo
    for idx = i_idx:f_idx
        # Resíduo da betweenness
        temp_array = abs.(betweenness[:, idx] - reference_betweenness)
        # Subtrai do mínimo
        temp_array = temp_array .- minimum(temp_array)
        # Normaliza para uma escala de 0 a 1
        temp_array = temp_array ./ maximum(temp_array)
        # Armazena a betwenness normalizada
        betweenness[:, idx] = temp_array
        for j = 1:n # Para cada nó
            # Se a variação da betweenness foi maior que 70%,
            # sinaliza como crítico
            if temp_array[j] > 0.7
                critical_nodes[j, idx] = 1
            elseif temp_array[j] < 0.2
                # Se a variação da betweenness foi menos que 20%,
                # sinaliza como não-critico
                noncritical_nodes[j, idx] = 1
            else
                medium_nodes[j, idx] = 1
            end
        end
    end
end

function myrange(q::SharedArray, dim)
    # Função auxiliar que distribui um SharedArray,
    # segundo sua dimensão dim, nos processos
    idx = indexpids(q)
    if idx == 0
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s)
              for s = range(0,
                            stop=size(q ,dim),
                            length=nchunks + 1)]
    splits[idx] + 1, splits[idx + 1]
end

s_current_flow_betweenness!(adjs::SharedArray,
                            betws::SharedArray
                            ) = current_flow_betweenness!(adjs::SharedArray,
                                                          betws::SharedArray,
                                                          myrange(adjs, 3))


s_verify_connectivity!(adjs::SharedArray,
                       tups::SharedArray
                       ) = verify_connectivity!(adjs::SharedArray,
                                                tups::SharedArray,
                                                myrange(adjs, 3))

s_edit_adjacencies!(adjs::SharedArray,
                    A_temp,
                    edg_tuples,
                    N::Int64
                    ) = edit_adjacencies!(adjs::SharedArray,
                                          A_temp,
                                          edg_tuples,
                                          N::Int64,
                                          myrange(adjs, 3))

s_define_critical_nodes!(nc_nodes::SharedArray,
                         md_nodes::SharedArray,
                         c_nodes::SharedArray,
                         betws::SharedArray,
                         ref_betws::SharedArray
                         ) = define_critical_nodes!(nc_nodes::SharedArray,
                                                    md_nodes::SharedArray,
                                                    c_nodes::SharedArray,
                                                    betws::SharedArray,
                                                    ref_betws::SharedArray,
                                                    myrange(c_nodes, 2))


"""
cfb_exaus(g, N)

Calcula a métrica para uma grafo `g`,
tomando remoções de `N` arestas por vez.
Foi definida pensando em redes grandes,
então possui parâmetros internos para limitar a
alocação de memória e realizar em loops o que
não for possível de ser realizado simultâneamente.
Sua implementação envolve o uso de funções
de computação paralela, mas a função não lida com a
criação de processos. Estes devem ser iniciados previamente.

### Parâmetros opcionais
- `chunk_limit=10000`: Limita a alocação máxima de memória das
matrizes de adjacência a uma variável de 10000MB. Valores maiores
demandam uma maior quantidade de RAM, mas agilizam o cálculo.
- `filename="current_flow"`: Parte do nome do arquivo `.txt`
que irá receber as informações finais.
"""
function cfb_exaus(g::Graph,
                   N::Int64,
                   chunk_limit=1000,
                   filename="result")
    println("Iniciando...")
     # Número de vértices
    n = nv(g)
    # Numero temporario de arestas
    m = ne(g)
    # ---- Variáveis importantes ----
    # Número total de remoções válidas
    total_removables = 0
    # Espaço de memória para a adjacência
    reference_adjacency = SharedArray{Bool}(n, n, 1)
    # Matriz de adjacência do grafo original
    reference_adjacency[:, :, 1] = adjacency_matrix(g)
    # Espaço de memória para a betweenness original
    reference_betweenness = SharedArray{Float64}(n, 1)
    # Cálculo da betweenness original
    current_flow_betweenness!(reference_adjacency,
                              reference_betweenness,
                              [1, 1])
    # Arestas do grafo original
    edgs = [[convert(UInt16, src(e)),
             convert(UInt16, dst(e))] for e = edges(g)]

    # Listas de 1 aresta 
    edgs_temp_singles =[e for e = subsets(edgs, 1)]
    # Memória compartilhada para as matrizes de adjacência
    adjacencies_temp = SharedArray{Bool}(n, n, m)
    tuples_to_remove_temp = SharedArray{Bool}(m)
    # Execução paralela da edição da matriz de adjacência para cada remoção
    println("Iniciando paralelismo")
    @sync begin
        for p = procs(adjacencies_temp)
            @async remotecall_wait(s_edit_adjacencies!,
                                   p,
                                   adjacencies_temp,
                                   reference_adjacency,
                                   edgs_temp_singles,
                                   1)
        end
    end
    # Execução paralela da verificação de conectividade para
    # cada matriz de adjacência
    println("Iniciando paralelismo")
    @sync begin
        for p = procs(tuples_to_remove_temp)
            @async remotecall_wait(s_verify_connectivity!,
                                   p,
                                   adjacencies_temp,
                                   tuples_to_remove_temp)
        end
    end
    # Considera apenas as matrizes dos grafos conexos
    valid_edgs = []
    for i in 1:m
        if tuples_to_remove_temp[i] == true
            append!(valid_edgs, [edgs[i]])
        end
    end
    # Coleta lixo para limpar RAM
    edgs_temp_singles = 0
    adjacencies_temp = 0
    tuples_to_remove_temp = 0
    GC.gc()
    # ---- Parâmetros utilizados na função ----
     # Número de arestas
    m = length(valid_edgs)
    # Número de grupos de N arestas existentes (binomial m N)
    tuples = multinomial(N, m - N)
    # Espaço de memória para armazenar tudo (MB)
    mem_chunk = 2 * 8 * n * n * tuples / 10^6
    println(string("Total de memória (MB): ", mem_chunk))
    # Número de iterações para executar tudo
    n_chunks = Int64(ceil(mem_chunk / chunk_limit))
    println(string("Número de chunks: ", n_chunks))
    # Dimensão dos vetores em cada iteração
    delta =  Int64(ceil(tuples / n_chunks))
    # Índices iniciais de cada iteração
    chunk_list = [1 + (i - 1) * delta for i = 1:n_chunks]
    # Listas de N arestas do grafo original
    edg_tuples = [e for e = subsets(valid_edgs, N)]
    # Listas com os índices das arestas de edg_tuples
    edg_tuples_indices = [convert(Array{UInt16}, i)
                          for i = subsets(1:m, N)]
    # Matrizes
    noncritical_matrix = SharedArray{Int32}(m, n)
    medium_matrix = SharedArray{Int32}(m, n)
    critical_matrix = SharedArray{Int32}(m, n)

    for chunk_iter = chunk_list
        println(string("Chunk ", chunk_iter))
        # Para armazenar as arestas cuja remoção foi válida
        edg_valid_removals = []
        # Verifica o limite superior da iteração atual
        ul = min(chunk_iter + delta - 1, tuples)
        # Calcula o delta real
        iter_delta = min(delta, ul - chunk_iter + 1)
        # Memória compartilhada para as matrizes de adjacência
        adjacencies = SharedArray{Bool}(n, n, iter_delta)
        # Vetor de bool para dizer quando uma remoção é válida (conexa)
        tuples_to_remove = SharedArray{Bool}(iter_delta)
        # Acessa as arestas necessárias para a iteração
        edg_tuples_chunk = edg_tuples[chunk_iter:ul]
        edg_tuples_chunk_indices = edg_tuples_indices[chunk_iter:ul]

        # Execução paralela da edição da matriz de adjacência para cada remoção
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(adjacencies)
                @async remotecall_wait(s_edit_adjacencies!,
                                       p,
                                       adjacencies,
                                       reference_adjacency,
                                       edg_tuples_chunk,
                                       N)
            end
        end
        # Execução paralela da verificação de conectividade
        # para cada matriz de adjacência
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(tuples_to_remove)
                @async remotecall_wait(s_verify_connectivity!,
                                       p,
                                       adjacencies,
                                       tuples_to_remove)
            end
        end
        # Contagem de quantas remoções foram válidas
        # (mantiveram o grafo conexo)
        removable = 0
        for i = 1:iter_delta
            if tuples_to_remove[i]
                removable += 1
                append!(edg_valid_removals,
                        [edg_tuples_chunk_indices[i]])
            end
        end
        total_removables += removable
        # Alocação de memória para as matrizes de adjacência
        # válidas e as betweenness
        final_adjacencies = SharedArray{Bool}(n, n, removable)
        betweenness = SharedArray{Float64}(n, removable)
        # Passagem das matrizes de adjacência válidas
        counter = 1
        for e = 1:iter_delta
            if tuples_to_remove[e]
                final_adjacencies[:, :, counter] = adjacencies[:, :, e]
                counter += 1
            end
        end
        # Execução paralela do cálculo das betweenness
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(betweenness)
                @async remotecall_wait(s_current_flow_betweenness!,
                                       p,
                                       final_adjacencies,
                                       betweenness)
            end
        end
        # Calcula a métrica e define os nós críticos
        noncritical_nodes = SharedArray{UInt16}(n, removable)
        medium_nodes = SharedArray{UInt16}(n, removable)
        critical_nodes = SharedArray{UInt16}(n, removable)
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(critical_nodes)
                @async remotecall_wait(s_define_critical_nodes!,
                                       p,
                                       noncritical_nodes,
                                       medium_nodes,
                                       critical_nodes,
                                       betweenness,
                                       reference_betweenness)
            end
        end
        # Atualiza a matriz
        for r = 1:removable
            for i = edg_valid_removals[r]
                for j = findall(noncritical_nodes[:, r] .!= 0)
                    noncritical_matrix[i, j] += 1
                end
                for j = findall(medium_nodes[:, r] .!= 0)
                    medium_matrix[i, j] += 1
                end
                for j = findall(critical_nodes[:, r] .!= 0)
                    critical_matrix[i, j] += 1
                end
            end
        end
		adjacencies = 0
		final_adjacencies = 0
		betweenness = 0
		noncritical_nodes = 0
		medium_nodes = 0
		critical_nodes = 0
        @everywhere GC.gc(true)
    end
    # Escreve os resultados em arquivos
    dir = string("exaustivo_", filename)
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    f = open("noncritical_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(noncritical_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("medium_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(medium_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("critical_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(critical_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("edges.txt", "w")
    for e in edgs
        write(f, string(e[1], ",", e[2], "\n"))
    end
    close(f)
    return "Sucesso!"
end


function cfb_de(g::Graph,
                N::Int64,
                pop_size::Int64,
                crossover_rate::Float64,
                iter_num=10::Int64)
    # Número de vértices
    n = nv(g)
    # Número de arestas
    m = ne(g)
    # ---- Variáveis importantes ----
    # Memória compartilhada para as matrizes de adjacência
    adjacencies = SharedArray{Int64}(n, n, pop_size)
    # Espaço de memória para a adjacência
    reference_adjacency = SharedArray{Bool}(n, n, 1)
    # Matriz de adjacência do grafo original
    reference_adjacency[:, :, 1] = adjacency_matrix(g)
    # Espaço de memória para a betweenness original 
    reference_betweenness = SharedArray{Float64}(n, 1)
    # Cálculo da betweenness original
    current_flow_betweenness!(reference_adjacency,
                              reference_betweenness,
                              [1, 1])
    # Arestas do grafo original
    edgs = [[src(e), dst(e)] for e = edges(g)]
    edgs_indices_pop = zeros(Int64, pop_size, N)
    # Parâmetros do DE
    beta_min = 0.2
    beta_max = 0.8
    # População inicial aleatória
    for i = 1:pop_size
        edgs_indices_pop[i, :] = sample(1:m, N, replace=false)    
    end
    best_cost = zeros(iter_num, 1)
    best_pop = []
    # Vetores de linhas críticas e nós críticos
    critical_lines = []
    critical_buses = []
    # Loop do algoritmo
    for iteration = 1:iter_num
        # Listas de N arestas do grafo original
        edg_tuples = [[edgs[edgs_indices_pop[i, j]]
                       for j in 1:N ] for i = 1:pop_size]
        # Vetor de bool para dizer quando uma remoção é válida (conexa)
        tuples_to_remove = SharedArray{Bool}(pop_size)
        # Execução paralela da edição da matriz de
        # adjacência para cada remoção
        @sync begin
            for p = procs(adjacencies)
                @async remotecall_wait(s_edit_adjacencies!,
                                       p,
                                       adjacencies,
                                       reference_adjacency,
                                       edg_tuples,
                                       N)
            end
        end
        # Execução paralela da verificação de conectividade
        # para cada matriz de adjacência
        @sync begin
            for p = procs(tuples_to_remove)
                @async remotecall_wait(s_verify_connectivity!,
                                       p,
                                       adjacencies,
                                       tuples_to_remove)
            end
        end
        # Contagem de quantas remoções foram válidas
        # (mantiveram o grafo conexo)
        removable = 0
        valid_edges_indices = []
        for i = 1:pop_size
            if tuples_to_remove[i]
                removable += 1
                append!(valid_edges_indices,
                        [edgs_indices_pop[i, :]])
            end
        end
        # Alocação de memória para as matrizes de
        # adjacência válidas e as betweenness
        final_adjacencies = SharedArray{Int64}(n, n, removable)
        betweenness = SharedArray{Float64}(n, removable)
        # Passagem das matrizes de adjacência válidas
        counter = 1
        for e = 1:pop_size
            if tuples_to_remove[e]
                final_adjacencies[:, :, counter] = adjacencies[:, :, e]
                counter += 1
            end
        end
        # Execução paralela do cálculo das betweenness
        @sync begin
            for p = procs(betweenness)
                @async remotecall_wait(s_current_flow_betweenness!,
                                       p,
                                       final_adjacencies,
                                       betweenness)
            end
        end
        # Calcula a métrica e define os nós críticos
        noncritical_nodes = SharedArray{UInt16}(n, removable)
        medium_nodes = SharedArray{UInt16}(n, removable)
        critical_nodes = SharedArray{Int64}(n, removable)
        @sync begin
            for p = procs(critical_nodes)
                @async remotecall_wait(s_define_critical_nodes!,
                                       p,
                                       noncritical_nodes,
                                       medium_nodes,
                                       critical_nodes,
                                       betweenness,
                                       reference_betweenness)
            end
        end
        # Obtém os custos e ordena os vetores de custos
        # e das arestas respectivas
        costs = reduce(+, critical_nodes, dims=[1])[1, :]
        p = sortperm(-costs)
        costs = costs[p]
        valid_edges_indices = valid_edges_indices[p]
        # Faz o backup dos melhores custos e as melhores
        # arestas associadas
        best_cost[iteration] = costs[1]
        append!(best_pop, [valid_edges_indices[1]])
        critical_buses_temp = findall(critical_nodes[:, p[1]] .!= 0)
        # Atualiza o conjunto de linhas críticas
        if iteration != 1
            if best_cost[iteration] > best_cost[iteration - 1]
                critical_lines = []
                critical_buses = []
            end
            i = 1
            while costs[i] == best_cost[iteration]
                critical_lines = union(critical_lines,
                                       [sort(valid_edges_indices[i])])
                critical_buses_temp = findall(critical_nodes[:, p[i]] .!= 0)
                critical_buses = union(critical_buses,
                                       [critical_buses_temp])
                i += 1
            end
        end

        # Realiza as operações de DE
        # Mutação
        coeffs = rand(-m:m, pop_size, N)
        beta = rand(pop_size, N).*(beta_max - beta_min) .+ beta_min
        edgs_indices_pop_mut = round.(Int64,
                                      edgs_indices_pop .+ (beta .* coeffs))
        edgs_indices_pop_mut = max.(ones(Int64, pop_size, N),
                                    edgs_indices_pop_mut)
        edgs_indices_pop_mut = min.(m*ones(Int64, pop_size, N),
                                    edgs_indices_pop_mut)
        # Crossover
        crossover_flag = rand(pop_size, 1) .< crossover_rate
        for c = 1:pop_size
            if crossover_flag[c]
                edgs_indices_pop[c, :] = edgs_indices_pop_mut[c, :]
            end
        end
        # Seleção simplificada - elitismo
        edgs_indices_pop[1, :] = valid_edges_indices[1]
    end
    println(string("Resultado: \n",
                   "     Melhor custo: ", best_cost[iter_num], "\n",
                   "     Linhas críticas: ", critical_lines, "\n",
                   "     Barras críticas: ", critical_buses)
           )
    return "Sucesso!"
end


"""
Calcula a mudança nas CFB em % da original, para
cada nó em cada remoção válida de tupla de arestas.
"""
function calculate_cfb_changes!(reference_betweenness::SharedArray,
                                betweenness::SharedArray,
                                betweenness_changes::SharedArray)
    println("Calculando mudanças na CFB")
    # Descobre o número de nós e de remoções válidas
    n, val_e = size(betweenness)
    for idx = 1:val_e
        # Resíduo da betweenness
        temp_array = abs.(betweenness[:, idx] - reference_betweenness)
        # Subtrai do mínimo
        temp_array = temp_array .- minimum(temp_array)
        # Normaliza para uma escala de 0 a 1
        temp_array = temp_array ./ maximum(temp_array)
        # Armazena a betwenness normalizada
        betweenness_changes[:, idx] = temp_array
    end
end


"""
cfb_guloso(g, N)

Calcula a métrica para uma grafo `g`,
tomando remoções de `N` arestas por vez.
Foi definida pensando em redes grandes,
então possui parâmetros internos para limitar a
alocação de memória e realizar em loops o que
não for possível de ser realizado simultâneamente.
Sua implementação envolve o uso de funções
de computação paralela, mas a função não lida com a
criação de processos. Estes devem ser iniciados previamente.

### Parâmetros opcionais
- `filename="current_flow"`: Parte do nome do arquivo `.txt`
que irá receber as informações finais.
"""
function cfb_guloso(g::Graph,
                    N::Int64,
                    filename="result")
    println("Iniciando...")
    # Número de vértices
    n = nv(g)
    # Numero temporario de arestas
    m = ne(g)
    # ---- Variáveis importantes ----
    # Número total de remoções válidas
    total_removables = 0
    # Espaço de memória para a adjacência
    reference_adjacency = SharedArray{Bool}(n, n, 1)
    # Matriz de adjacência do grafo original
    reference_adjacency[:, :, 1] = adjacency_matrix(g)
    # Espaço de memória para a betweenness original
    reference_betweenness = SharedArray{Float64}(n, 1)
    # Cálculo da betweenness original
    current_flow_betweenness!(reference_adjacency,
                              reference_betweenness,
                              [1, 1])
    # Arestas do grafo original
    edgs = [[convert(UInt16, src(e)),
             convert(UInt16, dst(e))] for e = edges(g)]

    # Listas de 1 aresta 
    edgs_temp_singles =[e for e = subsets(edgs, 1)]
    # Memória compartilhada para as matrizes de adjacência
    adjacencies_temp = SharedArray{Bool}(n, n, m)
    tuples_to_remove_temp = SharedArray{Bool}(m)
    # Execução paralela da edição da matriz de adjacência para cada remoção
    println("Iniciando paralelismo")
    @sync begin
        for p = procs(adjacencies_temp)
            @async remotecall_wait(s_edit_adjacencies!,
                                   p,
                                   adjacencies_temp,
                                   reference_adjacency,
                                   edgs_temp_singles,
                                   1)
        end
    end
    # Execução paralela da verificação de conectividade para
    # cada matriz de adjacência
    println("Iniciando paralelismo")
    @sync begin
        for p = procs(tuples_to_remove_temp)
            @async remotecall_wait(s_verify_connectivity!,
                                   p,
                                   adjacencies_temp,
                                   tuples_to_remove_temp)
        end
    end
    # Considera apenas as matrizes dos grafos conexos
    valid_edgs = []
    for i in 1:m
        if tuples_to_remove_temp[i] == true
            append!(valid_edgs, [edgs[i]])
        end
    end
    # Coleta lixo para limpar RAM
    edgs_temp_singles = 0
    adjacencies_temp = 0
    tuples_to_remove_temp = 0
    GC.gc()
    # ---- Parâmetros utilizados na função ----
     # Número de arestas
    m = length(valid_edgs)
    # Número de grupos de N arestas existentes (binomial m N)
    tuples = multinomial(N, m - N)
    # Listas de N arestas do grafo original
    edg_tuples = [e for e = subsets(valid_edgs, N)]
    # Listas com os índices das arestas de edg_tuples
    edg_tuples_indices = [convert(Array{UInt16}, i)
                          for i = subsets(1:m, N)]
    # Matrizes
    noncritical_matrix = SharedArray{Int32}(m, n)
    medium_matrix = SharedArray{Int32}(m, n)
    critical_matrix = SharedArray{Int32}(m, n)

    for chunk_iter = chunk_list
        println(string("Chunk ", chunk_iter))
        # Para armazenar as arestas cuja remoção foi válida
        edg_valid_removals = []
        # Verifica o limite superior da iteração atual
        ul = min(chunk_iter + delta - 1, tuples)
        # Calcula o delta real
        iter_delta = min(delta, ul - chunk_iter + 1)
        # Memória compartilhada para as matrizes de adjacência
        adjacencies = SharedArray{Bool}(n, n, iter_delta)
        # Vetor de bool para dizer quando uma remoção é válida (conexa)
        tuples_to_remove = SharedArray{Bool}(iter_delta)
        # Acessa as arestas necessárias para a iteração
        edg_tuples_chunk = edg_tuples[chunk_iter:ul]
        edg_tuples_chunk_indices = edg_tuples_indices[chunk_iter:ul]

        # Execução paralela da edição da matriz de adjacência para cada remoção
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(adjacencies)
                @async remotecall_wait(s_edit_adjacencies!,
                                       p,
                                       adjacencies,
                                       reference_adjacency,
                                       edg_tuples_chunk,
                                       N)
            end
        end
        # Execução paralela da verificação de conectividade
        # para cada matriz de adjacência
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(tuples_to_remove)
                @async remotecall_wait(s_verify_connectivity!,
                                       p,
                                       adjacencies,
                                       tuples_to_remove)
            end
        end
        # Contagem de quantas remoções foram válidas
        # (mantiveram o grafo conexo)
        removable = 0
        for i = 1:iter_delta
            if tuples_to_remove[i]
                removable += 1
                append!(edg_valid_removals,
                        [edg_tuples_chunk_indices[i]])
            end
        end
        total_removables += removable
        # Alocação de memória para as matrizes de adjacência
        # válidas e as betweenness
        final_adjacencies = SharedArray{Bool}(n, n, removable)
        betweenness = SharedArray{Float64}(n, removable)
        # Passagem das matrizes de adjacência válidas
        counter = 1
        for e = 1:iter_delta
            if tuples_to_remove[e]
                final_adjacencies[:, :, counter] = adjacencies[:, :, e]
                counter += 1
            end
        end
        # Execução paralela do cálculo das betweenness
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(betweenness)
                @async remotecall_wait(s_current_flow_betweenness!,
                                       p,
                                       final_adjacencies,
                                       betweenness)
            end
        end
        # Calcula a métrica e define os nós críticos
        noncritical_nodes = SharedArray{UInt16}(n, removable)
        medium_nodes = SharedArray{UInt16}(n, removable)
        critical_nodes = SharedArray{UInt16}(n, removable)
        println("Iniciando paralelismo")
        @sync begin
            for p = procs(critical_nodes)
                @async remotecall_wait(s_define_critical_nodes!,
                                       p,
                                       noncritical_nodes,
                                       medium_nodes,
                                       critical_nodes,
                                       betweenness,
                                       reference_betweenness)
            end
        end
        # Atualiza a matriz
        for r = 1:removable
            for i = edg_valid_removals[r]
                for j = findall(noncritical_nodes[:, r] .!= 0)
                    noncritical_matrix[i, j] += 1
                end
                for j = findall(medium_nodes[:, r] .!= 0)
                    medium_matrix[i, j] += 1
                end
                for j = findall(critical_nodes[:, r] .!= 0)
                    critical_matrix[i, j] += 1
                end
            end
        end
		adjacencies = 0
		final_adjacencies = 0
		betweenness = 0
		noncritical_nodes = 0
		medium_nodes = 0
		critical_nodes = 0
        @everywhere GC.gc(true)
    end

    # Escreve os resultados em arquivos
    dir = string("exaustivo_", filename)
    if !isdir(dir)
        mkdir(dir)
    end
    cd(dir)
    f = open("noncritical_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(noncritical_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("medium_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(medium_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("critical_matrix.txt", "w")
    for l = 1:m
        for c = 1:n
            write(f, string(critical_matrix[l, c], ","))
        end
        write(f, "\n")
    end
    close(f)
    f = open("edges.txt", "w")
    for e in edgs
        write(f, string(e[1], ",", e[2], "\n"))
    end
    close(f)
    return "Sucesso!"
end


"""
iteracao_cfb(g)

Realiza uma procura pela aresta única cuja remoção promove a
maior diferença em termos de CFB num grafo  `g`.

Retorna a aresta cuja remoção foi a mais crítica.
"""
function iteracao_cfb(g::Graph)
    println("Iniciando...")
    # Número de vértices
    n = nv(g)
    # Numero temporario de arestas
    m = ne(g)
    # ---- Variáveis importantes ----
    # Número total de remoções válidas
    total_removables = 0
    # Espaço de memória para a adjacência
    reference_adjacency = SharedArray{Bool}(n, n, 1)
    # Matriz de adjacência do grafo original
    reference_adjacency[:, :, 1] = adjacency_matrix(g)
    # Espaço de memória para a betweenness original
    reference_betweenness = SharedArray{Float64}(n, 1)
    # Cálculo da betweenness original
    current_flow_betweenness!(reference_adjacency,
                              reference_betweenness,
                              [1, 1])
    
    # Gera a lista de arestas do grafo
    edgs = [[convert(UInt16, src(e)),
             convert(UInt16, dst(e))] for e = edges(g)]
    edgs_temp_singles = [e for e = subsets(edgs, 1)]
    calc_range = [1, length(edgs_temp_singles)]
    # Gera as matrizes de adjacências editadas
    adjacencies_temp = SharedArray{Bool}(n, n, m)
    edit_adjacencies!(adjacencies_temp,
                      reference_adjacency,
                      edgs_temp_singles,
                      1,
                      calc_range)

    # Verifica as matrizes de grafos conexos
    tuples_to_remove_temp = SharedArray{Bool}(m)
    verify_connectivity!(adjacencies_temp,
                         tuples_to_remove_temp,
                         calc_range)

    # Contagem de quantas remoções foram válidas
    # (mantiveram o grafo conexo)
    valid_edgs = []
    valid_edgs_indices = []
    for i = 1:m
        if tuples_to_remove_temp[i] == true
            append!(valid_edgs, [edgs[i]])
            append!(valid_edgs_indices, [i])
        end
    end

    removable = length(valid_edgs)
    final_adjacencies = SharedArray{Bool}(n, n, removable)
    betweenness = SharedArray{Float64}(n, removable)
    # Passagem das matrizes de adjacência válidas
    counter = 1
    for e = 1:removable
        if tuples_to_remove_temp[e]
            final_adjacencies[:, :, counter] = adjacencies_temp[:, :, e]
            counter += 1
        end
    end
    # Calcula as betweenness para cada remoção
    current_flow_betweenness!(final_adjacencies,
                              betweenness,
                              [1, removable])
    betweenness_changes = SharedArray{Float64}(n, removable)
    calculate_cfb_changes!(reference_betweenness,
                           betweenness,
                           betweenness_changes)
    # Obtém a remoção mais crítica (com mais nós críticos)
    most_critical = 0
    most_critical_edge_idx = 0
    for e = 1:removable
        critical_removal = length(findall(betweenness_changes[:, e] .> 0.7))
        if critical_removal > most_critical
            most_critical = critical_removal
            most_critical_edge_idx  = e
        end
    end
    return valid_edgs[most_critical_edge_idx], most_critical, betweenness_changes
end

end
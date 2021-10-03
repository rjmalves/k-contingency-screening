import networkx as nx
from networkx.algorithms.components.connected import is_connected
import numpy as np
from copy import deepcopy
from itertools import combinations

arq = "t11.txt"
k = 1

grafo = nx.read_edgelist(arq)
n = grafo.number_of_nodes()
m = grafo.number_of_edges()
cfb_referencia = nx.current_flow_betweenness_centrality(grafo)

# Monta as tuplas de arestas
arestas = list(grafo.edges)
tuplas = list(combinations(arestas, k))

# Gera os grafos com tuplas removidas
grafos_editados = []
for t in tuplas:
    grafo_copia = deepcopy(grafo)
    for aresta in t:
        grafo_copia.remove_edge(*aresta)
    grafos_editados.append(grafo_copia)


# Verifica remoções de tuplas válidas
grafos_conexos = []
tuplas_validas = []
for g, t in zip(grafos_editados, tuplas):
    if nx.is_connected(g):
        grafos_conexos.append(g)
        tuplas_validas.append(t)

n_tuplas_validas = len(tuplas_validas)

# Monta a matriz de deltas para as tuplas válidas
matriz_deltas = np.zeros((n, n_tuplas_validas))

# Daqui pra baixo tem que fazer (:
# Preenche a matriz de deltas
for i, g in enumerate(grafos_conexos):
    # Calcula cfb do grafo com a tupla removida
    cfb_grafo = nx.current_flow_betweenness_centrality(g)
    # Preenche na t-ésima coluna da matriz
    matriz_deltas[:, i] = list(cfb_grafo.values())

# Resume em um valor global por aresta

# k-contingency-screening
Code for analyzing multiple contingencies in Power Systems by graph theory in Julia.


## Documentação dos Resultados

Os arquivos de saída fornecidos pelo método exaustivo e pela evolução diferencial (DE) serão descritos a seguir.

### disconnects.csv

Contém as tuplas de arestas que desconectam o grafo quando removidas. Para estar, não é possível quantificar o impacto, pois ele é "infinito".
Possui `2k` colunas, onde `k` é o número de arestas em cada tupla. As colunas (1, 2), (3, 4), e assim sucessivamente representam os vértices de origem e destino das arestas de cada tupla.

### valid_tuples.csv

Contém as tuplas de arestas que não desconectam o grafo e, portanto, são ditas "válidas". Assim como o arquivo anterior, possui `2k` colunas, onde `k` é o número de arestas de cada tupla e cada par de colunas contém os vértices de origem e destino das arestas da tupla.

### local_deltas.csv

Contém o impacto local da remoção de cada tupla de arestas nos vértices. Este impacto é calculado através do somatório das diferenças da métrica em questão, para todos os vértices, numa mesma remoção. Estes somatórios (um para cada tupla), são então normalizados para o intervalo [0, 1]. O valor da linha `i` contém o impacto local da remoção da tupla da linha `i` do arquivo `valid_tuples.csv`.

### vertex_global_deltas.csv

Contém o impacto global nos vértices causado pela remoção das tuplas. Este impacto é calculado através da soma das diferenças da métrica em questão (somatório dos deltas) para um mesmo vértice, após considerar a remoção de todas as tuplas, uma a uma. Isto é chamado de `impacto global nos vértices`. Os valores escritos no arquivo são os desta soma, normalizada para o intervalo de 0 a 1. A linha `i` contém o valor do impacto global para o vértice `i`.

### edge_global_deltas.csv

Contém o impacto global nas arestas, causado pela remoção das tuplas. Este impacto é calculado, para cada aresta `e`, através da soma dos impactos locais (antes de normalizá-los), presentes no arquivo `local_deltas.csv` em todas as tuplas nas quais a aresta `e` está presente. A linha `i` contém o impacto global da remoção da aresta definida na linha `i` no arquivo de lista de arestas do grafo.



# k-contingency-screening

Código para análise topológica de contingências múltiplas em grafos simples que modelam sistemas de potência, utilizando a Current-Flow Betweenness Centrality.


## Configuração e Visão Geral

O repositório consiste em em três variantes principais de métodos para realizar a mesma análise: análise e ordenação e contingências modeladas como remoção de arestas em grafos.

São implementados os métodos:

- Exaustivo
- Guloso
- Evolução Diferencial (DE)

As funcionalidades são ilustradas no arquivo `script.jl`, que pode ser livremente modificado.

Uma vez clonado o repositório na máquina local e possuindo uma instalação de Julia (recomendado > 1.9), a configuração do ambiente pode ser feita através de:

```
$ cd <caminho_do_repositorio>
$ julia
$ ]
$ activate .
$ instantiate
```

O comando `]` é utilizado para entrar no modo de gerenciamento de pacotes de Julia. Mais informações podem ser encontradas em https://docs.julialang.org/en/v1/stdlib/Pkg/ .

Dentro da própria sessão interativa de julia, estando no ambiente em que foram instaladas as dependências (basta conferir o que aparece entre parênteses na linha do terminal quando se digita `]`), o conteúdo do script pode ser executado com:

```
$ include("script.jl")
```

Para executar diretamente do terminal, fora de uma sessão interativa de Julia, basta fazer, estando dentro do diretório em que foi clonado o repositório:

```
julia --project script.jl
```

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

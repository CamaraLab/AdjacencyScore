## Installation

```
install.packages("expm")
devtools::install_github("CamaraLab/AdjacencyScore", auth_token = _________ )
library(AdjacencyScore)
```

## Functions

```
knn_graph(data, k=5)
```
Calculates an approximate knn graph using RANN and outputs the adjacency matrix in sparse matrix format.

```
adjacency_score(adj_matrix, f, f_pairs, k)
```
Returns the adjacency score, p value, and q value for each pair of features in f_pairs.

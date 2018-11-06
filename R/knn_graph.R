#' Returns the adjacency matrix of a knn graph from the given distance matrix
#'
#' @param data a matrix of points as rows
#' @param k number of nearest neighbors to compute
#'
#' @import RANN
#' @import igraph
#' @export
#'


knn_graph <- function(data, k) {
  knn <- nn2(data, k=k)$nn.idx
  el <- cbind(knn[, 1], c(knn[, -1]))
  adj <- get.adjacency(graph.edgelist(el))
  sym_adj <- 1*((adj + t(adj)) > 0)
  return(sym_adj)
}

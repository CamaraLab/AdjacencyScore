% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlap_score.R
\name{overlap_score}
\alias{overlap_score}
\title{Ranks pairs of features by how much they overlap in a graph}
\usage{
overlap_score(f, f_pairs, num_perms = 1000, seed = 10, num_cores = 1,
  perm_estimate = F)
}
\arguments{
\item{f}{a numeric matrix specifying one or more features defined for each node of the graph.
Each column is a node of the graph and each row is a feature over the nodes.}

\item{f_pairs}{a 2 column matrix where each row specifies the indices or names
of a pair of features on which the score will be computed}

\item{num_perms}{number of permutations used to build the null distribution for each
feature. By default is set to 1000.}

\item{seed}{integer specifying the seed used to initialize the generator of permutations.
By default is set to 10.}

\item{num_cores}{integer specifying the number of cores to be used in the computation. By
default only one core is used.}

\item{perm_estimate}{boolean indicating whether normal distribution parameters should be
determined from num_perms permutations to estimate the p-value. By default is set to FALSE.}
}
\description{
Given the adjacency matrix of a graph and a set of features on that graph, ranks given pairs
of those features (f and g) by the equation fg, which measures how much those features
cooccur in the same cells. Calculates the p-value for this score by permuting
the columns of the feature matrix separately for each features.
}

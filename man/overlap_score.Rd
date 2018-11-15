% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlap_score.R
\name{overlap_score}
\alias{overlap_score}
\title{Ranks pairs of features by how much they overlap in a graph}
\usage{
overlap_score(f, f_pairs, num_perms = 1000, seed = 10, num_cores = 1)
}
\arguments{
\item{f}{a numeric vector or matrix specifying one or more functions with support on
the set of points whose significance will be assesed in the simplicial complex. Each
column corresponds to a point and each row specifies a different function.}

\item{f_pairs}{a 2 column matrix where each row specifes the indices or names
of a pair of points on which the Comb. Lap. score will be computed}

\item{num_perms}{number of permutations used to build the null distribution for each
feature. By default is set to 1000.}

\item{seed}{integer specifying the seed used to initialize the generator of permutations.
By default is set to 10.}

\item{num_cores}{integer specifying the number of cores to be used in the computation. By
default only one core is used.}
}
\description{
Given the adjacency matrix of a graph and a set of features on that graph, ranks given pairs
of those features (f and g) by the equation fg, which measures how much those features
cooccur in the same cells. Calculates the p-value for this score by permuting
the columns of the feature matrix separately for each features.
}
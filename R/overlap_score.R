#' Ranks pairs of features using the Combinatorial Laplacian Score for 0- and 1-forms.
#'
#' Given a nerve or a clique complex, a set of features consisting of functions with support on
#' the set of points underlying the complex, and a list of pairs of features,
#' it asseses the significance of each pair of features
#' in the simplicial complex by computing its scalar and vectorial Combinatorial Laplacian
#' Score and comparing it with the null distribution that results from reshufling many times the values of
#' the function across the point cloud. For nerve complexes, feature functions induce 0- and
#' 1-forms in the complex by averaging the function across the points associated to 0- and 1-simplices
#' respectively. For clique complexes, feature functions are directly 0-forms in the complex and 1-forms
#' are obtained by averaging the function across the two vertices connected by each edge.
#'
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param f_pairs a 2 column matrix where each row specifes the indices or names
#' of a pair of points on which the Comb. Lap. score will be computed
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#'
#' @import Matrix
#' @import parallel
#' @import expm
#' @import MASS
#' @export

overlap_score <- function(f, f_pairs, num_perms = 1000, seed = 10, num_cores = 1) {

  # Check class of f
  if (class(f) != 'matrix') {
    f <- as.matrix(f)
  }

  # Check class of f_pairs
  if (class(f_pairs) == 'list') {
    f_pairs <- matrix(unlist(f_pairs), ncol=2, byrow=T)
  } else if (class(f_pairs) == 'numeric' || class(f_pairs) == 'character') {
    # only one pair
    f_pairs <- matrix(f_pairs, ncol=2, byrow=T)
  }

  qqh <- NULL
  for (i in 1:nrow(f_pairs)) {
    f1 <- as.numeric(f[f_pairs[i,1],])
    f2 <- as.numeric(f[f_pairs[i,2],])
    if (num_perms > 0) {
      perm1 <- t(mcmapply(function(x) f1[sample(1:ncol(f))], 1:num_perms, mc.cores=num_cores))
      perm2 <- t(mcmapply(function(x) f2[sample(1:ncol(f))], 1:num_perms, mc.cores=num_cores))
      f1 <- rbind(f1, perm1)
      f2 <- rbind(f2, perm2)
      qt <- rowSums(f1*f2)
    } else {
      qt <- sum(f1 * f2)
    }
    ph <- c( score = qt[1], p = (sum(qt>=qt[1])-1.0)/num_perms)
    qqh <- rbind(qqh, ph)
  }
  qqh <- as.data.frame(qqh, row.names=1:nrow(qqh))

  # Adjust for multiple hypothesis testing
  qqh$q <- p.adjust(qqh$p, method = 'BH')

  # Add feature columns
  qqh <- data.frame(f = f_pairs[,1], g = f_pairs[,2], qqh, stringsAsFactors=F)

  return(qqh)
}

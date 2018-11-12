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
#' @param g2 an object of the class \code{simplicial} containing the nerve or clique complex.
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param f_pairs a 2 column matrix where each row specifes the indices or names
#' of a pair of points on which the Comb. Lap. score will be computed
#' @param c constant used to determine width of diffusion, must be 0 <= c
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
#' @export

adjacency_score <- function(adj_matrix, f, f_pairs, c, num_perms = 1000, seed = 10, num_cores = 1) {

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

  permutations <- t(mcmapply(function(x) sample(1:ncol(f)), 1:num_perms, mc.cores=num_cores))
  permutations <- rbind(1:ncol(f), permutations)

  # Permute and normalize each feature, result is a list of matrices where each matrix corresponds to all the permutations for each feature
  perm_f <- mclapply(1:nrow(f), function(i) t(sapply(1:nrow(permutations), function(j) f[i,][permutations[j,]] - sum(f[i,])/ncol(f))), mc.cores=num_cores)
  names(perm_f) <- row.names(f)

  adj_sym <- 1*((adj_matrix+t(adj_matrix)) > 0)
  diag(adj_sym) <- 0

  if (c != 0) {
    expm_adj <- expm::expm(c*adj_sym, method="Higham08")
    expm_adj <- as.matrix(expm_adj)
    expm_adj <- (expm_adj - diag(nrow(adj_sym)))/c
  }

  # Evaluates R and p for a pair of features fo
  cornel <- function(fo) {
    f1 <- perm_f[[fo[1]]]
    f2 <- perm_f[[fo[2]]]
    if (c == 0) {
      qt <- rowSums((f1%*%adj_sym)*f2)
    } else {
      qt <- rowSums((f1%*%expm_adj)*f2)
    }
    ph <- NULL
    ph$score <- qt[1]
    ph$p <- (sum(qt>=qt[1])-1.0)/num_perms
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of pairs of features
  worker <- function(fu) {
    qh <- NULL
    qh$score <- NULL
    qh$p <- NULL
    for (i in 1:nrow(fu)) {
      d <- cornel(fu[i,])
      qh$score <- rbind(qh$score, d$score)
      qh$p <- rbind(qh$p, d$p)
    }
    return(data.frame(qh))
  }

  if (num_cores > nrow(f_pairs)) {
    num_cores <- nrow(f_pairs)
  }

  if (num_cores == 1 || nrow(f_pairs) == 1) {
    qqh <- worker(f_pairs)
  } else {
    # If more than one core then split the features in num_cores parts accordingly
    wv <- floor(nrow(f_pairs)/num_cores)
    wr <- nrow(f_pairs) - wv*num_cores
    work <- list()
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (f_pairs[(1+(m-1)*(wv+1)):(m*(wv+1)),,drop=FALSE])
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (f_pairs[(1+wr+(m-1)*wv):(wr+m*wv),,drop=FALSE])
      }
    } else {
      for (m in 1:num_cores) {
        work[[m]] <- (f_pairs[(1+(m-1)*wv):(m*wv),,drop=FALSE])
      }
    }
    reul <- mclapply(work, worker, mc.cores = num_cores)
    qqh <- reul[[1]]
    for (m in 2:num_cores) {
      qqh <- rbind(qqh, reul[[m]])
    }
  }

  # Adjust for multiple hypothesis testing
  qqh$q <- p.adjust(qqh$p, method = 'BH')

  # Add feature columns
  qqh <- data.frame(f = f_pairs[,1], g = f_pairs[,2], qqh, stringsAsFactors=F)

  return(qqh)
}

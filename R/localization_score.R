#' Ranks features by how localized they are in graph
#'
#' Given the adjacency matrix of a graph and a set of features on that graph, ranks those features (f)
#' by the equation f((e^{kA}-I)/k)f, which measures how much those features are localized in the graph.
#' Calculates the p-value of this score by permuting the columns of the feature matrix.
#'
#' @param adj_matrix a (preferrably sparse) binary matrix of adjacency between the columns of f
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param perm_estimate boolean indicating whether normal distribution parameters should be
#' determined from num_perms permutations to estimate the p-value. By default is set to FALSE.
#'
#' @import Matrix
#' @import parallel
#' @import MASS
#' @export

localization_score <- function(adj_matrix, f, num_perms = 1000, num_cores = 1, perm_estimate = F) {

  # Check class of f
  if (class(f) != 'matrix') {
    f <- as.matrix(f)
  }

  permutations <- NULL
  if (num_perms > 0) {
    permutations <- t(mcmapply(function(x) sample(1:ncol(f)), 1:num_perms, mc.cores=num_cores))
  }

  permutations <- rbind(1:ncol(f), permutations)

  # Permute and normalize each feature, result is a list of matrices where each matrix corresponds to all the permutations for each feature
  perm_f <- mclapply(1:nrow(f), function(i) t(sapply(1:nrow(permutations), function(j) f[i,][permutations[j,]])), mc.cores=num_cores)
  names(perm_f) <- row.names(f)

  adj_sym <- 1*((adj_matrix+t(adj_matrix)) > 0)
  diag(adj_sym) <- 0

  # Evaluates R and p for feature fo
  cornel <- function(fo) {
    qt <- rowSums((fo%*%adj_sym)*fo)
    ph <- NULL
    ph$score <- qt[1]
    if (perm_estimate) {
      nfit <- fitdistr(as.numeric(qt), "normal")
      ph$p <- 1-pnorm(qt[1], mean=nfit$estimate["mean"], sd = nfit$estimate["sd"])
    }
    else {
      ph$p <- sum(qt>qt[1])/num_perms
    }
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of features
  worker <- function(fu) {
    qh <- NULL
    qh$score <- NULL
    qh$p <- NULL
    for (i in 1:length(fu)) {
      d <- cornel(fu[[i]])
      qh$score <- rbind(qh$score, d$score)
      qh$p <- rbind(qh$p, d$p)
    }
    return(data.frame(qh))
  }

  if (num_cores > length(perm_f)) {
    num_cores <- length(perm_f)
  }

  if (num_cores == 1 || length(perm_f) == 1) {
    qqh <- worker(perm_f)
  } else {
    # If more than one core then split the features in num_cores parts accordingly
    wv <- floor(length(perm_f)/num_cores)
    wr <- length(perm_f) - wv*num_cores
    work <- list()
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (perm_f[(1+(m-1)*(wv+1)):(m*(wv+1))])
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (perm_f[(1+wr+(m-1)*wv):(wr+m*wv)])
      }
    } else {
      for (m in 1:num_cores) {
        work[[m]] <- (perm_f[(1+(m-1)*wv):(m*wv)])
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

  row.names(qqh) <- row.names(f)

  return(qqh)
}

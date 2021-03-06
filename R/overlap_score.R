#' Ranks pairs of features by how much they overlap in a graph
#'
#' Given the adjacency matrix of a graph and a set of features on that graph, ranks given pairs
#' of those features (f and g) by the equation fg, which measures how much those features
#' cooccur in the same cells. Calculates the p-value for this score by permuting
#' the columns of the feature matrix separately for each features.
#'
#' @param f a numeric matrix specifying one or more features defined for each node of the graph.
#' Each column is a node of the graph and each row is a feature over the nodes.
#' @param f_pairs a 2 column matrix where each row specifies the indices or names
#' of a pair of features on which the score will be computed
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param perm_estimate boolean indicating whether normal distribution parameters should be
#' determined from num_perms permutations to estimate the p-value. By default is set to FALSE.
#'
#' @import Matrix
#' @import parallel
#' @import MASS
#' @export

overlap_score <- function(f, f_pairs, num_perms = 1000, seed = 10, num_cores = 1, perm_estimate = F) {

  # Check class of f
  if (!is(f,'matrix') && !is(f,'Matrix')) {
    cat("Converting f to matrix\n")
    f <- as.matrix(f)
  }
  
  # Check class of f_pairs
  if (is(f_pairs,'list')) {
    f_pairs <- matrix(unlist(f_pairs), ncol=2, byrow=T)
  } else if (is(f_pairs,'numeric') || is(f_pairs,'character')) {
    # only one pair
    f_pairs <- matrix(f_pairs, ncol=2, byrow=T)
  }

  perm_f <- vector('list', nrow(f))
  for (i in 1:nrow(f)) {
    permutations <- NULL
    if (num_perms > 0) {
      permutations <- t(mcmapply(function(x) sample(1:ncol(f)), 1:num_perms, mc.cores=num_cores))
    }
    permutations <- rbind(1:ncol(f), permutations)
    perm_f[[i]] <-  t(sapply(1:nrow(permutations), function(j) f[i,][permutations[j,]]))
  }
  names(perm_f) <- row.names(f)

  # Evaluates R and p for a pair of features fo
  cornel <- function(fo) {
    f1 <- perm_f[[fo[1]]]
    f2 <- perm_f[[fo[2]]]
    if (num_perms > 0) {
      qt <- rowSums(f1*f2)
    } else {
      qt <- sum(f1 * f2)
    }
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

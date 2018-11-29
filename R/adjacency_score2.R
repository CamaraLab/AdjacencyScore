#' Ranks pairs of features by how localized they are in graph
#'
#' Given the adjacency matrix of a graph and a set of features on that graph, ranks given pairs
#' of those features (f and g) by the equation f((e^{kA}-I)/k)g, which measures how much those
#' features are colocalized in the graph. Calculates the p-value for this score by permuting
#' the columns of the feature matrix separately for each features.
#'
#' @param adj_matrix a (preferrably sparse) binary matrix of adjacency between the columns of f
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param f_pairs a 2 column matrix where each row specifes the indices or names
#' of a pair of points on which the Comb. Lap. score will be computed
#' @param c constant used to determine width of diffusion, must be 0 <= c
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param perm_estimate boolean indicating whether normal distribution parameters should be
#' determined from num_perms permutations to estimate the p-value. By default is set to FALSE.
#' @param groupings boolean indicating whether features are binary and mutually exclusive
#' indicated each point's inclusion in some group. Allows for p-value computation from a
#' parameterized hypergeometric null distribution. By default is set to FALSE.
#'
#' @import Matrix
#' @import parallel
#' @import expm
#' @import MASS
#' @export

adjacency_score2 <- function(adj_matrix, f, f_pairs, c, num_perms = 1000, num_cores = 1, perm_estimate = F, groupings=F) {

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

  if (c != 0 && groupings) {
    groupings <- FALSE
    cat("Setting groupings to FALSE since c > 0\n")
  }

  # If using parameterized distribution for p-value, don't need permutations
  if (groupings && num_perms > 0) {
    num_perms <- 0
    cat("Setting num_perms to 0 since using grouping features\n")
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
    if (groupings) {
      # features are mutually exclusive
      overlap <- sum(f1 * f2)
      if (overlap == 0) {
        wb <- 2 * sum(f1) * sum(f2)
        edges <- sum(adj_sym)/2
        ph$p <- 1-phyper(qt[1], m=wb, n = length(f1)*(length(f1)-1) - wb, k = edges)
      }
      # features entirely overlap (same grouping)
      else {
        wb <- overlap*(overlap - 1)
        edges <- sum(adj_sym)/2
        ph$p <- 1-phyper(floor(qt[1]/2), m=wb, n = length(f1)*(length(f1)-1) - wb, k = edges)
      }
    }
    else if (perm_estimate) {
      nfit <- fitdistr(as.numeric(qt), "normal")
      ph$p <- 1-pnorm(qt[1], mean=nfit$estimate["mean"], sd = nfit$estimate["sd"])
    }
    else {
      ph$p <- (sum(qt>=qt[1])-1.0)/num_perms
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

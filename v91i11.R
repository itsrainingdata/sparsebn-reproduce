### Package installation ---------------------------------------------
# Make sure that the package and its dependencies have been installed
pkgs <- c("sparsebn", "discretecdAlgorithm", "ccdrAlgorithm", "sparsebnUtils", "network", "BiocManager", "SID")

install_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}
instl <- sapply(pkgs, install_if_needed)

if (!requireNamespace("Rgraphviz", quietly = TRUE)) BiocManager::install("Rgraphviz")

### Options ----------------------------------------------------------
options(digits = 3)  # to suppress output of excessive decimal places
set.seed(1)  # for reproducibility

### Section 1 --------------------------------------------------------
library("sparsebn")
data("pathfinder")
data <- sparsebnData(pathfinder$data, type = "continuous")
dags <- estimate.dag(data)
plotDAG(dags)

# Re-run with timer
system.time({
  dags <- estimate.dag(data)
})

# Figure 1: Plot first 4 non-trivial DAGs in solution path
# Note: First estimate is always the null graph
plotDAG(dags[2:5])

### Section 4 --------------------------------------------------------

### Section 4.1 --------------------------------------------------------
library("sparsebn")
library("bnlearn")
library("Matrix")
library("mvtnorm")
library("pcalg")

gen_block_dag <- function(nblocks = 2) {
  ### Load DAG
  load(url("http://www.bnlearn.com/bnrepository/pathfinder/pathfinder.rda"))
  gr <- as.graphNEL(bn)
  pathfinderDAG <- to_edgeList(as.graphNEL(bn))

  ### Construct block diagonal matrix from DAG
  coefs <- get.adjacency.matrix(pathfinderDAG)
  # nblocks <- 5
  coefs <- Matrix::bdiag(rep(list(coefs), nblocks))

  coefs
}

gen_gaussian_data <- function(coefs, nn = 50) {
  nnode <- ncol(coefs)

  ### Randomly flip half the signs
  nedge <- sum(coefs > 0)
  prop.flip <- 0.5
  flipped <- sample(which(coefs > 0), floor(prop.flip * nedge))
  coefs[flipped] <- -1 * coefs[flipped]

  ### Generate random data
  vars <- Matrix::Diagonal(n = nnode, x = rep(1, nnode))
  id <- vars
  covMat <- t(solve(id - coefs)) %*% vars %*% solve(id - coefs)
  gaussian.data <- rmvnorm(n = nn, sigma = as.matrix(covMat))

  sparsebnData(gaussian.data, type = "c")
}

run_scalability_test <- function(nblocks, nn = 50, algs = c("ccdr", "pc", "mmhc")) {
  stopifnot(is.list(nblocks))

  timing <- vector("list", length = length(nblocks))
  for (k in seq_along(nblocks)) {
    dat <- gen_gaussian_data(coefs = gen_block_dag(nblocks[[k]]), nn = nn)

    times <- rep(NA, 8)
    edges <- rep(NA, 4)
    names(edges) <- c("ccdr", "pc", "mmhc", "hc")

    ### Comparison times
    if ("pc" %in% algs) {
      library("pcalg")
      times[3] <- proc.time()[3]
      out <- pc(suffStat = list(C = cor(dat$data), n = nrow(dat$data)), indepTest = gaussCItest,
        alpha = 0.01, labels = names(dat$data), verbose = FALSE)
      times[4] <- proc.time()[3]

      edges["pc"] <- graph::numEdges(out@graph)
    }

    if ("mmhc" %in% algs) {
      library("bnlearn")
      times[5] <- proc.time()[3]
      out <- mmhc(dat$data)
      times[6] <- proc.time()[3]

      edges["mmhc"] <- nrow(out$arcs)
    }

    if ("hc" %in% algs) {
      times[7] <- proc.time()[3]
      out <- hc(dat$data)
      times[8] <- proc.time()[3]

      edges["hc"] <- nrow(out$arcs)
    }

    # Set edge threshold for CCDr
    threshold <- NA
    threshold <- max(edges, na.rm = T)
    if (is.na(threshold))
      threshold <- 5000

    if ("ccdr" %in% algs) {
      library("sparsebn")
      times[1] <- proc.time()[3]
      dags.out <- estimate.dag(data = dat, lambdas.length = 20, edge.threshold = threshold,
        verbose = FALSE)
      times[2] <- proc.time()[3]

      edges["ccdr"] <- max(num.edges(dags.out))
    }

    timing[[k]] <- list(nblocks = nblocks[[k]], pp = ncol(dat$data), nn = nrow(dat$data),
      ccdr = times[2] - times[1], pc = times[4] - times[3], mmhc = times[6] -
        times[5], hc = times[8] - times[7], threshold = threshold, ccdr = edges["ccdr"],
      pc = edges["pc"], mmhc = edges["mmhc"], hc = edges["hc"])
    print(do.call("rbind", timing))
  }

  do.call("rbind", timing)
}

### Run simulation (takes ~10 minutes)
pathfinder <- run_scalability_test(as.list(1:10), 50)

### Reproduce Figure 2
dat <- pathfinder[, c(2, 4, 5, 6)]
matplot(dat[, 1], dat[, -1], pch = c("C", "P", "M"), xlab = "p (number of nodes)",
  ylab = "Time (s)", type = "b", col = c("black", "blue", "red"))

### Section 4.2 --------------------------------------------------------
library("igraph")

########################################################################
# useful functions
########################################################################

##################
# get_summary       a function to get summary for every dag along solution path
# adjMatrix.path    a list of every adjacency matrix along the solution path
# adjMatrix.true    the true adjacency matrix
# rdm_order         a random order to permutate data set
# return            a list of summaries for every graph along solution path
get_summary <- function(adjMatrix.path, adjMatrix.true, rdm_order) {
  sapply(adjMatrix.path, function(adjMatrix, adjMatrix.true, rdm_order) {
    temp_adjMatrix <- matrix(0, ncol = ncol(adjMatrix.true), nrow = nrow(adjMatrix.true))
    temp_adjMatrix[rdm_order, rdm_order] <- as.matrix(adjMatrix)
    adjMatrix <- temp_adjMatrix
    # get P, E, R, FP
    P <- sum(adjMatrix)
    E <- sum(adjMatrix * adjMatrix.true)
    R <- sum(t(adjMatrix) * adjMatrix.true)
    FP <- P - E - R
    # get SHD
    shd <- SID::hammingDist(adjMatrix, adjMatrix.true)
    return(c(P, E, R, FP, shd))
  }, adjMatrix.true, rdm_order)
}

##################
# run_cd              a function that run and time cd algorithm or adaptive cd algorithm
# databn              a sparsebnData object
# data_type data type that needed to be generated, 'continuous' or 'discrete'
# adjMatrix.true      adjacency matrix of true graph
# adapt               TRUE or FALSE, if TRUE both adaptive lasso and lasso will be run.
# rdm_order           a random order to permutate data set
# lambda.seq          lambda sequence for regular cd algorithm
# lambda.ratio        lambda ratio for regular cd algorithm
# lambda.length       length of lambda sequence for regular cd algorithm
# lambda.scale        method to scale lambda sequence for regular cd algorithm. Has 'linear' and 'log'
# ...                 See estimate.dag function and discretecdAlgorithm::cd.run function
# return              solution path for both cd algorithm or adaptive cd algorithm.
#                     brief summary for every dag along the solution path.
run_cd <- function(databn, data_type, adjMatrix.true, adapt = FALSE, rdm_order, ...) {
  if (data_type == "continuous") {
    begin <- Sys.time()
    cd.path <- sparsebn::estimate.dag(databn, ...)
    end <- Sys.time()
  } else if (data_type == "discrete") {
    begin <- Sys.time()
    cd.path <- discretecdAlgorithm::cd.run(databn, adaptive = adapt, ...)
    end <- Sys.time()
  } else {
    stop("Please specify a valid data type: continuous or discrete!!!")
  }

  # get running time of cd algorithm
  # calculate time in seconds
  time <- as.numeric(end - begin, units = "secs")

  # get adjacency matrix for CD algorithm solution path
  matrix.cd.path <- sparsebnUtils::get.adjacency.matrix(cd.path)

  # get summary for cd algorithm
  cd.summary <- get_summary(matrix.cd.path, adjMatrix.true, rdm_order)
  cd.summary <- as.data.frame(t(cd.summary))
  cd.summary <- cbind(cd.summary, rep(time, nrow(cd.summary)))
  colnames(cd.summary) <- c("P", "E", "R", "FP", "SHD", "time")
  rownames(cd.summary) <- NULL

  return(list(cd.path = cd.path, summary = cd.summary, time = time))
}

##################
# match_edges     a function that return index of a dag that has a number of edges that's closest to n_edge
# cd.path         a sparsebnPath object
# n_edge          target number of edges
# return          an index of a dag that has closest number of edges compared with n_edge
match_edges <- function(cd.path, n_edge) {
  cd_nedges <- sapply(cd.path, function(x) {
    x$nedge
  })
  index <- which.min(abs(cd_nedges - n_edge))
  if (index == 1) {
    if (length(cd.path) >= 2) {
      index <- 2
    }
  }

  return(index)
}

##################
# data_prep       a function to prepare data set
# data.intervene  data set that has one intervention on a variable for each observation
# data.observe    observational data set
# n_ivn           number of interventions needed
# n_obs           number of total observations
# node_index      a vector of nodes should be under interventions
# return          a sparsebnData object
data_prep <- function(data.intervene, data.observe, n_ivn, n_obs, node_index = NULL) {
  if ((ncol(data.intervene) - 1) != ncol(data.observe))
    stop("Data sets not compatible!! \n Check input data!")
  if (n_ivn * length(node_index) > n_obs)
    stop("n_ivn is too big!")
  if (n_obs <= 0)
    stop("Observation should be a positive number!")
  if (n_ivn < 0)
    stop("number of interventions should be a non-negative number!")

  n_node <- ncol(data.intervene) - 1

  ivn_total <- data.intervene[, ncol(data.intervene)]
  data_itv <- data.intervene[, -ncol(data.intervene)]
  data <- c()

  if (is.null(node_index)) {
    data <- data.observe[1:n_obs, ]
    # ivn <- lapply(1:n_obs, function(x){0L})
    ivn <- lapply(1:n_obs, function(x) {
      NULL
    })
    databn <- sparsebnUtils::sparsebnData(data, ivn = ivn, type = DATA_TYPE)

    return(databn)
  }

  if (n_ivn > 0) {
    for (i in 1:length(node_index)) {
      single_index <- node_index[i]
      data <- rbind(data, data_itv[which(ivn_total == single_index), ][1:n_ivn,
        ])
    }
    if ((n_obs - length(node_index) * n_ivn) > 0) {
      data <- rbind(data, data.observe[1:(n_obs - length(node_index) * n_ivn),
        ])
    }
  } else {
    data <- data.observe[1:n_obs, ]
  }

  ivn_vec <- c(rep(node_index, rep(n_ivn, length(node_index))), rep(0, (n_obs -
    length(node_index) * n_ivn)))
  ivn <- lapply(ivn_vec, function(x) {
    as.integer(x)
  })
  for (i in 1:n_obs) {
    if (ivn[[i]] == 0) {
      ivn[i] <- list(NULL)
    }
  }


  databn <- sparsebnUtils::sparsebnData(data, ivn = ivn, type = DATA_TYPE)
  return(databn)
}

##################
# cd_intervene_all    function that test the influence of
#                     increading number of interventions on each node.
# data.intervene      interventional data set. Last column should be index for intervention.
# data.observe        observational data set
# adjMatrix.true      adjacency matrix for true graph
# k_max               maximum number of interventions on each node
# n_obs               total number of interventions needed
# rdm_order           a random order to permutate data set
# verbose             TRUE or FALSE value
# return              a list of summary tables for cd algorithm and adaptove cd algorithm.
cd_intervene_all <- function(data.intervene, data.observe, adjMatrix.true, k_max,
  n_obs, rdm_order, verbose = TRUE) {
  n_node <- ncol(data.intervene) - 1
  T <- sum(adjMatrix.true)
  node_index <- 1:n_node

  if (n_node * k_max > n_obs) {
    stop("number of intervention is too big. k_max * number of nodes must be smaller or equals to n_obs!")
  }

  all.est <- vector("list", 1)
  names(all.est) <- c("CD")
  # names(all.est) <- c("CD", "CDAdaptive")
  for (i in 1:(k_max + 1)) {
    k <- i - 1
    if (verbose == TRUE)
      cat("Running data with number of interventions:", k, "... \n")
    databn <- data_prep(data.intervene, data.observe, k, n_obs, node_index)
    cd.out <- run_cd(databn, data_type = DATA_TYPE, adjMatrix.true = adjMatrix.true,
      rdm_order = rdm_order)
    # cdAdapt.out <- run_cd(databn, data_type = DATA_TYPE, adjMatrix.true = adjMatrix.true, adapt = TRUE)
    index <- match_edges(cd.out$cd.path, T)
    # indexAdapt <- match_edges(cdAdapt.out$cd.path, T)
    temp_summary <- cd.out$summary[index, ]
    row.names(temp_summary) <- NULL
    # tempAdapt_summary <- cdAdapt.out$summary[indexAdapt, ]
    # row.names(tempAdapt_summary) <- NULL
    all.est$CD <- rbind(all.est$CD, temp_summary)
    # all.est$CDAdaptive <- rbind(all.est$CDAdaptive, tempAdapt_summary)
  }
  all.est <- lapply(all.est, as.data.frame)

  summary.all <- lapply(all.est, function(x, n_node, T) {
    k <- 0:k_max
    P <- x$P
    E <- x$E
    R <- x$R
    FP <- x$FP
    SHD <- x$SHD
    time <- x$time
    TPR <- x$E/T
    FDR_P <- x$P
    FDR_P[which(x$P == 0)] <- 1
    FDR <- (x$R + x$FP)/FDR_P
    FPR <- (x$R + x$FP)/(1/2 * n_node * (n_node - 1) - T)

    summary <- cbind(P = P, E = E, R = R, FP = FP, SHD = SHD, TPR = TPR, FDR = FDR,
      FPR = FPR, time = time)
    row.names(summary) <- k
    return(summary)
  }, n_node, T)

  return(summary.all)
}

##################
# cd_multi_test       function that test only part of the nodes are under intervention
# data.intervene      interventional data set. Last column should be index for intervention.
# data.observe        observational data set
# adjMatrix.true      adjacency matrix for true graph
# n_ivn               number of interventions for each node
# n_obs               total number of interventions needed
# rdm_order           a random order to permutate data set
# return              a list of summary tables for cd algorithm and adaptove cd algorithm.
cd_multi_test <- function(databn, adjMatrix.true, n_ivn, n_obs, rdm_order) {
  n_node <- ncol(databn$data)
  T <- sum(adjMatrix.true)

  # all.est <- vector("list", 2)
  # names(all.est) <- c("CD", "CDAdaptive")
  all.est <- vector("list", 1)
  names(all.est) <- c("CD")

  # node_index <- sort(node_index)
  cd.out <- run_cd(databn, data_type = DATA_TYPE, adjMatrix.true = adjMatrix.true,
    rdm_order = rdm_order)
  # cdAdapt.out <- run_cd(databn, data_type = DATA_TYPE, adjMatrix.true = adjMatrix.true, adapt = TRUE)
  index <- match_edges(cd.out$cd.path, T)
  # indexAdapt <- match_edges(cdAdapt.out$cd.path, T)
  temp_summary <- cd.out$summary[index, ]
  row.names(temp_summary) <- NULL
  # tempAdapt_summary <- cdAdapt.out$summary[indexAdapt, ]
  # row.names(tempAdapt_summary) <- NULL
  all.est$CD <- rbind(all.est$CD, temp_summary)
  # all.est$CDAdaptive <- rbind(all.est$CDAdaptive, tempAdapt_summary)
  all.est <- lapply(all.est, as.data.frame)

  summary.all <- lapply(all.est, function(x, n_node, T) {
    P <- x$P
    E <- x$E
    R <- x$R
    FP <- x$FP
    SHD <- x$SHD
    time <- x$time
    TPR <- x$E/T
    FDR_P <- x$P
    FDR_P[which(x$P == 0)] <- 1
    FDR <- (x$R + x$FP)/FDR_P
    FPR <- (x$R + x$FP)/(1/2 * n_node * (n_node - 1) - T)

    summary <- cbind(P = P, E = E, R = R, FP = FP, SHD = SHD, TPR = TPR, FDR = FDR,
      FPR = FPR, time = time)
    return(summary)
  }, n_node, T)

  return(summary.all)
}

##################
# cd_loop             function that test the influence of
#                     increasing number of nodes under intervention
# data.intervene      interventional data set. Last column should be index for intervention.
# data.observe        observational data set
# adjMatrix.true      adjacency matrix for true graph
# n_ivn               number of interventions for each node
# n_obs               total number of interventions needed
# rdm_order           a random order to permutate data set
# track               TRUE or FALSE to keep track of number of tests
# return              a list of summary tables for cd algorithm and adaptove cd algorithm.
cd_loop <- function(data.intervene, data.observe, adjMatrix.true, n_ivn, n_obs, rdm_order,
  track = FALSE) {
  cd_total <- c()
  cdAdapt_total <- c()

  n_node <- ncol(data.intervene) - 1
  max_ivn <- n_node

  for (i in 0:max_ivn) {
    if (track == TRUE)
      cat("running tests with number of nodes under intervention is:", i, "...\n")
    n <- i
    all_nodes <- 1:n_node
    if (n == 0) {
      node_index <- NULL
    } else if (is.null(node_index)) {
      node_index <- sample(all_nodes, 1)
    } else {
      if (length(all_nodes[-node_index]) > 1) {
        node_index <- c(node_index, sample(all_nodes[-node_index], 1))
      } else {
        node_index <- c(node_index, all_nodes[-node_index])
      }
    }

    databn <- data_prep(data.intervene, data.observe, n_ivn, n_obs, node_index)

    summary <- cd_multi_test(databn, adjMatrix.true, n_ivn, n_obs, rdm_order = rdm_order)
    cd_total <- rbind(cd_total, summary$CD)
    # cdAdapt_total <- rbind(cdAdapt_total, summary$CDAdaptive)
  }
  rownames(cd_total) <- (0:max_ivn)
  # rownames(cdAdapt_total) <- (0:max_ivn)

  cd_total <- as.data.frame(cd_total)
  # cdAdapt_total <- as.data.frame(cdAdapt_total)

  # return(list(cd_total = cd_total, cdAdapt_total = cdAdapt_total))
  return(list(cd_total = cd_total))
}

##################
# cd_sub_run          function that test the influence of
#                     increasing number of nodes under interventions, average of n_tests
# data.intervene      interventional data set. Last column should be index for intervention.
# data.observe        observational data set
# adjMatrix.true      adjacency matrix for true graph
# n_ivn               number of interventions for each node
# n_obs               total number of interventions needed
# n_tests             number of tests
# rdm_order           a random order to permutate data set
# verbose             TRUE or FALSE value to keep track of the number of sets of tests.
# track               TRUE or FALSE value to keep track of tests in each set of tests
# return              a list of summary tables for cd algorithm and adaptove cd algorithm.
cd_sub_run <- function(data.intervene, data.observe, adjMatrix.true, n_ivn, n_obs,
  n_tests, rdm_order, verbose = TRUE, track = FALSE) {
  n_node <- ncol(data.intervene) - 1
  if (n_ivn * n_node > n_obs) {
    stop("n_ivn too big, n_ivn * number of nodes must be smaller or equals to n_obs!")
  }

  cd_average <- matrix(0, (n_node + 1), 9)
  cd_average <- as.data.frame(cd_average)
  # cdAdapt_average <- matrix(0, (n_node+1), 9); cdAdapt_average <- as.data.frame(cdAdapt_average)

  for (i in 1:n_tests) {
    if (verbose == TRUE) {
      cat("running the", i, "th set of tests \n")
    }
    summary <- cd_loop(data.intervene = data.intervene, data.observe = data.observe,
      adjMatrix.true = adjMatrix.true, n_ivn = n_ivn, n_obs = n_obs, rdm_order = rdm_order,
      track = track)
    cd_average <- cd_average + summary$cd_total
    # cdAdapt_average <- cdAdapt_average+summary$cdAdapt_total
  }
  cd_average <- cd_average/n_tests
  cd_average <- as.data.frame(cd_average)
  # cdAdapt_average <- cdAdapt_average/n_tests; cdAdapt_average <- as.data.frame(cdAdapt_average)
  colnames(cd_average) <- c("P", "E", "R", "FP", "SHD", "TPR", "FDR", "FPR", "time")
  # colnames(cdAdapt_average) <- c("P", "E", "R", "FP", "SHD", "TPR", "FDR", "FPR", "time")

  # return(list(cd_average = cd_average, cdAdapt_average = cdAdapt_average))
  return(list(cd_average = cd_average))
}

##################
# test_coef_gen  generate coefficients for binary data
#
# graph          igraph object
# return         coefficient list, assuming every variable is binary
test_coef_gen <- function(graph) {
  edge_list <- as.list(sparsebnUtils::to_edgeList(graph))
  coef <- vector("list", length = length(edge_list))
  coef <- lapply(edge_list, function(x) {
    sub_coef <- NULL
    if (length(x)) {
      sub_coef <- matrix(c(2, 0, rep(c(-2, 2), length(x))), nrow = 2)
    }
    sub_coef
  })
}

##################
# gen_sim_data Method to simulate data for experiments
# graph        igraph object
# n.obs        number of observational data set
# n.ivn        number of interventional data set
# ivn          List of intervention
# rdm_order    a random order to permutate data set
# type         Type of data set, 'discrete' or 'continuous'
gen_sim_data <- function(graph, n.obs, n.ivn, ivn, rdm_order, type = "discrete") {
  ivn_name_list <- lapply(ivn, function(x) {
    paste0("V", x)
  })
  if (type == "discrete") {
    el <- sparsebnUtils::as.edgeList(graph)
    names(el) <- paste0("V", 1:sparsebnUtils::num.nodes(el))
    coef <- test_coef_gen(graph)

    # observational data
    data.obs <- discretecdAlgorithm::generate_discrete_data(graph = el, n = n,
      params = coef)

    # interventional data
    data.ivn <- discretecdAlgorithm::generate_discrete_data(graph = el, n = n,
      params = coef, ivn = ivn_name_list)
    # data.ivn <- cbind(data.ivn, rep(1:p, rep(n/p, p)))
    # data.ivn <- cbind(data.ivn, unlist(ivn)) # last column is the list of intervention
  } else if (type == "continuous") {
    el <- sparsebnUtils::as.edgeList(graph)
    names(el) <- paste0("V", 1:sparsebnUtils::num.nodes(el))
    coef <- sparsebnUtils:::gen_params(el)

    # ivn_name_list <- lapply(ivn, function(x){paste0("V", x)})
    # observational data
    data.obs <- ccdrAlgorithm::generate_mvn_data(graph = el, n = n, params = coef)

    # interventional data
    data.ivn <- ccdrAlgorithm::generate_mvn_data(graph = el, n = n, params = coef,
      ivn = ivn_name_list)
    # data.ivn <- cbind(data.ivn, rep(1:p, rep(n/p, p))) # last column is the list of intervention
    # data.ivn <- cbind(data.ivn, unlist(ivn)) # last column is the list of intervention
  } else {
    stop()
  }

  data.obs <- as.data.frame(data.obs)
  data.ivn <- as.data.frame(data.ivn)
  data.obs <- data.obs[, rdm_order]
  data.ivn <- data.ivn[, rdm_order]
  rdm_ivn <- lapply(ivn_name_list, function(x) match(x, colnames(data.ivn)))
  data.ivn <- cbind(data.ivn, ivn = unlist(rdm_ivn))
  list(obs = data.obs, ivn = data.ivn)
}

########################################################################
# graph generation
########################################################################

set.seed(1)

### Set this flag to specify which type of data to test
DATA_TYPE <- "discrete"
# DATA_TYPE <- "continuous"

# set parameters
p <- 50  # number of parameters
n <- 500  # number of observations

#####################
# generate bipartite graph
p1 <- 10  # number of edges only have in edges
p2 <- 40  # number of edges only have out edges
graph_bipartite <- sample_bipartite(n1 = p1, n2 = p2, m = p, type = "gnm", directed = TRUE)

# generate polytree
n_subtree <- 5
p_subtree <- p/n_subtree
graph_polytree <- make_tree(n = p_subtree, mode = "out")
# build up unconnected trees
for (i in 2:n_subtree) {
  graph_polytree <- graph_polytree %du% make_tree(n = p_subtree, mode = "out")
}
# add edge between sub-trees
graph_polytree <- add_edges(graph_polytree, c(8, 31))
graph_polytree <- add_edges(graph_polytree, c(19, 32))
graph_polytree <- add_edges(graph_polytree, c(20, 41))
graph_polytree <- add_edges(graph_polytree, c(27, 41))

# generate scale-free graph
graph_scalefree <- sample_pa(n = p, m = 1, directed = FALSE)
edge_scalefree <- as_edgelist(graph_scalefree)
graph_scalefree <- graph(rbind(edge_scalefree[, 1], edge_scalefree[, 2]))

# generate small-world graph
graph_smallworld <- sample_smallworld(dim = 1, size = p, nei = 2, p = 0.05)
edge_smallworld <- as_edgelist(graph_smallworld)
graph_smallworld <- graph(rbind(edge_smallworld[, 1], edge_smallworld[, 2]))

########################################################################
# data generation
########################################################################

# list of intervention
ivn_list <- as.list(rep(1:p, rep(n/p, p)))
rdm_order <- sample(1:p, p)
# ivn_name_list <- lapply(ivn_list, function(x){paste0("V", x)})

data.bipartite <- gen_sim_data(graph = graph_bipartite, n.obs = n, n.ivn = n, ivn = ivn_list,
  rdm_order = rdm_order, type = DATA_TYPE)
data.polytree <- gen_sim_data(graph = graph_polytree, n.obs = n, n.ivn = n, ivn = ivn_list,
  rdm_order = rdm_order, type = DATA_TYPE)
data.scalefree <- gen_sim_data(graph = graph_scalefree, n.obs = n, n.ivn = n, ivn = ivn_list,
  rdm_order = rdm_order, type = DATA_TYPE)
data.smallworld <- gen_sim_data(graph = graph_smallworld, n.obs = n, n.ivn = n, ivn = ivn_list,
  rdm_order = rdm_order, type = DATA_TYPE)

########################################################################
# true DAG adjacency matrix generation
########################################################################

adjMatrix.true.bipartite <- as.matrix(sparsebnUtils::get.adjacency.matrix(sparsebnUtils::to_edgeList(graph_bipartite)))
adjMatrix.true.polytree <- as.matrix(sparsebnUtils::get.adjacency.matrix(sparsebnUtils::to_edgeList(graph_polytree)))
adjMatrix.true.scalefree <- as.matrix(sparsebnUtils::get.adjacency.matrix(sparsebnUtils::to_edgeList(graph_scalefree)))
adjMatrix.true.smallworld <- as.matrix(sparsebnUtils::get.adjacency.matrix(sparsebnUtils::to_edgeList(graph_smallworld)))

########################################################################
# FIGURE 3: Test the influence of increasing number of interventions on each node
########################################################################
k_max <- 10
n_obs <- 500

summary_all.scalefree <- cd_intervene_all(data.scalefree$ivn, data.scalefree$obs,
  adjMatrix.true.scalefree, k_max, n_obs, rdm_order = rdm_order, verbose = TRUE)
summary_all.smallworld <- cd_intervene_all(data.smallworld$ivn, data.smallworld$obs,
  adjMatrix.true.smallworld, k_max, n_obs, rdm_order = rdm_order, verbose = TRUE)
summary_all.polytree <- cd_intervene_all(data.polytree$ivn, data.polytree$obs, adjMatrix.true.polytree,
  k_max, n_obs, rdm_order = rdm_order, verbose = TRUE)
summary_all.bipartite <- cd_intervene_all(data.bipartite$ivn, data.bipartite$obs,
  adjMatrix.true.bipartite, k_max, n_obs, rdm_order = rdm_order, verbose = TRUE)
sum_all_cd_sf <- as.data.frame(summary_all.scalefree$CD)
sum_all_cd_sw <- as.data.frame(summary_all.smallworld$CD)
sum_all_cd_pt <- as.data.frame(summary_all.polytree$CD)
sum_all_cd_bp <- as.data.frame(summary_all.bipartite$CD)

# plot
par(cex.lab = 1.5)
x <- 0:10
TPR_low <- c(min(c(sum_all_cd_sf$TPR, sum_all_cd_sw$TPR, sum_all_cd_pt$TPR, sum_all_cd_bp$TPR)))
TPR_high <- 1

# a plot for Figure 3A
par(mfrow = c(1, 1), oma = c(4, 4, 1, 1) + 0.5, mar = c(1, 1, 1, 1) + 0.8)
plot(c(0, 10), c(TPR_low, TPR_high), pch = " ")
title(xlab = "m (number of interventions per node)", ylab = "TPR", outer = TRUE,
  line = 1)
lines(x, sum_all_cd_sf$TPR, lty = 1, col = "green")
lines(x, sum_all_cd_sw$TPR, lty = 5, col = "red")
lines(x, sum_all_cd_pt$TPR, lty = 6, col = "blue")
lines(x, sum_all_cd_bp$TPR, lty = 3, col = "black")

# a plot for Figure 3B
par(mfrow = c(2, 2), oma = c(5, 4, 0, 0) + 0.5, mar = c(1, 1, 1, 1) + 0.8)
plot(x, sum_all_cd_sf$TPR, pch = " ", main = "Scale-free")
lines(x, sum_all_cd_sf$TPR, lty = 1, col = "green")
plot(x, sum_all_cd_sw$TPR, pch = " ", main = "Small-world")
lines(x, sum_all_cd_sw$TPR, lty = 5, col = "red")
plot(x, sum_all_cd_pt$TPR, pch = " ", main = "Polytree")
lines(x, sum_all_cd_pt$TPR, lty = 6, col = "blue")
plot(x, sum_all_cd_bp$TPR, pch = " ", main = "Bipartite")
lines(x, sum_all_cd_bp$TPR, lty = 3, col = "black")
title(xlab = "m (number of interventions per node)", ylab = "TPR", outer = TRUE,
  line = 2)

################################################################################
# FIGURE 4: Test the influence of increasing number of nodes under intervention
################################################################################

n_ivn <- 10
n_obs <- 500
n_tests <- 20  # Number of sets of tests for one graph. Each set will have 50 individual tests.
# When setting n_tests to be 20, like what we did in the paper,
# please expect the running time to be more than a few hours.

summary.scalefree <- cd_sub_run(data.scalefree$ivn, data.scalefree$obs, adjMatrix.true.scalefree,
  n_ivn, n_obs, n_tests, rdm_order = rdm_order, verbose = TRUE, track = TRUE)
summary.smallworld <- cd_sub_run(data.smallworld$ivn, data.smallworld$obs, adjMatrix.true.smallworld,
  n_ivn, n_obs, n_tests, rdm_order = rdm_order, verbose = TRUE, track = TRUE)
summary.polytree <- cd_sub_run(data.polytree$ivn, data.polytree$obs, adjMatrix.true.polytree,
  n_ivn, n_obs, n_tests, rdm_order = rdm_order, verbose = TRUE, track = TRUE)
summary.bipartite <- cd_sub_run(data.bipartite$ivn, data.bipartite$obs, adjMatrix.true.bipartite,
  n_ivn, n_obs, n_tests, rdm_order = rdm_order, verbose = TRUE, track = TRUE)
sum_cd_sf <- summary.scalefree$cd_average
sum_cd_sw <- summary.smallworld$cd_average
sum_cd_pt <- summary.polytree$cd_average
sum_cd_bp <- summary.bipartite$cd_average

# plot
x <- 0:50
TPR_low <- c(min(c(sum_cd_sf$TPR, sum_cd_sw$TPR, sum_cd_pt$TPR, sum_cd_bp$TPR)))
TPR_high <- 1

# a plot for Figure 4A
par(mfrow = c(1, 1), oma = c(4, 4, 1, 1) + 0.5, mar = c(1, 1, 1, 1) + 0.8)
plot(c(0, 50), c(TPR_low, TPR_high), pch = " ")
title(xlab = "k (number of nodes under intervention)", ylab = "TPR", outer = TRUE,
  line = 1)
lines(x, sum_cd_sf$TPR, lty = 1, col = "green")
lines(x, sum_cd_sw$TPR, lty = 5, col = "red")
lines(x, sum_cd_pt$TPR, lty = 6, col = "blue")
lines(x, sum_cd_bp$TPR, lty = 3, col = "black")

# a plot for Figure 4B
par(mfrow = c(2, 2), oma = c(5, 4, 0, 0) + 0.5, mar = c(1, 1, 1, 1) + 0.8)
plot(x, sum_cd_sf$TPR, pch = " ", xlab = "k", ylab = "TPR", main = "Scale-free")
lines(x, sum_cd_sf$TPR, lty = 1, col = "green")
plot(x, sum_cd_sw$TPR, pch = " ", xlab = "k", ylab = "TPR", main = "Small-world")
lines(x, sum_cd_sw$TPR, lty = 5, col = "red")
plot(x, sum_cd_pt$TPR, pch = " ", xlab = "k", ylab = "TPR", main = "Polytree")
lines(x, sum_cd_pt$TPR, lty = 6, col = "blue")
plot(x, sum_cd_bp$TPR, pch = " ", xlab = "k", ylab = "TPR", main = "Bipartite")
lines(x, sum_cd_bp$TPR, lty = 3, col = "black")
title(xlab = "m (number of interventions per node)", ylab = "TPR", outer = TRUE,
  line = 2)

### Section 4.6 --------------------------------------------------------
# install.packages("sparsebn")
# devtools::install_github(c("itsrainingdata/sparsebnUtils@dev",
#                            "itsrainingdata/ccdrAlgorithm@dev",
#                            "gujyjean/discretecdAlgorithm@dev",
#                            "itsrainingdata/sparsebn@dev"))

### Section 5 --------------------------------------------------------
library("sparsebn")
data("cytometryContinuous")
names("cytometryContinuous")

### Section 5.1 --------------------------------------------------------
cyto.raw <- cytometryContinuous$data
cyto.ivn <- cytometryContinuous$ivn
cyto.data <- sparsebnData(cyto.raw, type = "continuous", ivn = cyto.ivn)
print(cyto.data)
summary(cyto.data)

### Section 5.2 --------------------------------------------------------
cyto.learn <- estimate.dag(cyto.data)
print(cyto.learn)
summary(cyto.learn)

estimate.dag(cyto.data, lambdas.length = 50)

### Use a linear scale
cyto.lambdas <- generate.lambdas(lambda.max = 10, lambdas.ratio = 0.001, lambdas.length = 10,
  scale = "linear")
cyto.lambdas

### Use a log scale
cyto.lambdas <- generate.lambdas(lambda.max = 10, lambdas.ratio = 0.001, lambdas.length = 10,
  scale = "log")
cyto.lambdas

estimate.dag(cyto.data, lambdas = cyto.lambdas)

### Section 5.3 --------------------------------------------------------
whitelist <- matrix(c("pip3", "pip2"), nrow = 1)
estimate.dag(cyto.data, whitelist = whitelist)

whitelist

blacklist <- rbind(c("raf", "jnk"), c("jnk", "raf"))
estimate.dag(cyto.data, blacklist = blacklist)

blacklist <- specify.prior(roots = "pip3", leaves = "jnk", nodes = names(cyto.data$data))
estimate.dag(cyto.data, blacklist = blacklist)

### Section 5.4 --------------------------------------------------------
print(cyto.learn[[1]])
summary(cyto.learn[[1]])

print(cyto.learn[[3]])
summary(cyto.learn[[3]])

get.adjacency.matrix(cyto.learn[[3]])

show.parents(cyto.learn[[3]], c("raf", "pip2"))

### Section 5.5 --------------------------------------------------------
cyto.param <- estimate.parameters(cyto.learn, data = cyto.data)

cyto.param[[3]]$coefs

Matrix::diag(cyto.param[[3]]$vars)

### Section 5.6 --------------------------------------------------------
select(cyto.learn, edges = 8)  # exact match returned
select(cyto.learn, edges = 10)  # closest match returned

select(cyto.learn, lambda = 41.75)  # closest match returned
select(cyto.learn, lambda = 41.7)  # same output as previous line

select(cyto.learn, index = 4)  # exact match returned
cyto.learn[[4]]  # same output as previous line

selected.lambda <- select.parameter(cyto.learn, cyto.data)
selected.lambda

### Section 5.7 --------------------------------------------------------
getPlotPackage()

setPlotPackage("network")
getPlotPackage()

setPlotPackage("graph")
getPlotPackage()

setPlotPackage("igraph")
getPlotPackage()

plot(cyto.learn)

plot(cytometryContinuous[["dag"]], layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
  igraph::in_circle()), vertex.label = names(cytometryContinuous[["dag"]]), vertex.size = 30,
  vertex.label.color = gray(0), vertex.color = gray(0.9), edge.color = gray(0),
  edge.arrow.size = 0.5)
plot(cyto.learn[[selected.lambda]], layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
  igraph::in_circle()), vertex.label = get.nodes(cyto.learn), vertex.size = 30,
  vertex.label.color = gray(0), vertex.color = gray(0.9), edge.color = gray(0),
  edge.arrow.size = 0.5)

# Figure 7 (code not in paper)
data(cytometryContinuous)
plot(cytometryContinuous[["dag"]],
     layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
                              igraph::in_circle()),
     vertex.label = names(cytometryContinuous[["dag"]]),
     vertex.size = 30,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.5
)

cyto.raw <- cytometryContinuous[["data"]]
cyto.ivn <- cytometryContinuous[["ivn"]]
cyto.data <- sparsebnData(cytometryContinuous[["data"]],
                          type = "continuous",
                          ivn = cytometryContinuous[["ivn"]])
cyto.learn <- estimate.dag(cyto.data)
plot(cyto.learn[[8]],
     layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
                              igraph::in_circle()),
     vertex.label = get.nodes(cyto.learn),
     vertex.size = 30,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.5
)

data(cytometryDiscrete)
cyto.data <- sparsebnData(cytometryDiscrete[["data"]],
                                 type = "discrete",
                                 ivn = cytometryDiscrete[["ivn"]])
cyto.learn <- estimate.dag(data = cyto.data)
plot(select(cyto.learn, edges = 19),
     layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
                              igraph::in_circle()),
     vertex.label = get.nodes(cyto.learn),
     vertex.size = 30,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.5
)

# Figure 6 (code not in paper)
par(mfrow = c(1, 3))
setPlotPackage("igraph")
plot(cyto.learn[[selected.lambda]])
setPlotPackage("network")
plot(cyto.learn[[selected.lambda]])
setPlotPackage("graph")
plot(cyto.learn[[selected.lambda]])

setPlotPackage("igraph")  # reset to default

### The code below requires the RCy3 package to be installed
### along with the Cytoscape application (http://www.cytoscape.org/).
###
### Cytoscape must be installed separately.
if (!requireNamespace("RCy3", quietly = TRUE)) BiocManager::install("RCy3")

openCytoscape(cytometryContinuous$dag)
openCytoscape(cyto.learn[[selected.lambda]])

### Section 6.1 --------------------------------------------------------
library("sparsebn")
data("cytometryDiscrete")

cyto.data <- sparsebnData(cytometryDiscrete$data, type = "discrete", ivn = cytometryDiscrete$ivn)
cyto.data

cyto.learn <- estimate.dag(cyto.data)
cyto.learn

plot(select(cyto.learn, edges = 19), layout = igraph::layout_(to_igraph(cytometryContinuous$dag),
  igraph::in_circle()), vertex.label = get.nodes(cyto.learn), vertex.size = 30,
  vertex.label.color = gray(0), vertex.color = gray(0.9), edge.color = gray(0),
  edge.arrow.size = 0.5)

cyto.param <- estimate.parameters(cyto.learn, data = cyto.data)
cyto.param[[5]][["raf"]]

### Section 6.2 --------------------------------------------------------
data("pathfinder")
dat <- sparsebnData(pathfinder$data, type = "c")

nn <- num.samples(dat)
lambdas <- generate.lambdas(sqrt(nn), 0.05, lambdas.length = 50, scale = "linear")
dags <- estimate.dag(data = dat, lambdas = lambdas, edge.threshold = 500, verbose = FALSE)
dags

dat <- sparsebnData(pathfinder$data[1:50, ], type = "c")
nn <- num.samples(dat)
lambdas <- generate.lambdas(sqrt(nn), 0.05, lambdas.length = 50, scale = "linear")
dags <- estimate.dag(data = dat, lambdas = lambdas, edge.threshold = 500, verbose = FALSE)
dags

### Section 6.3 --------------------------------------------------------
### The same code snippet can be run on all three datasets below
library("sparsebn")

path <- "./alzheimers44768PFC_5000.csv"
geneexpr <- read.csv(path, check.names = FALSE)
dat <- sparsebnData(geneexpr, type = "c")
lambdas <- generate.lambdas(sqrt(nrow(geneexpr)), scale = "linear", lambdas.length = 20)
system.time({
  dags <- estimate.dag(dat, lambdas = lambdas, edge.threshold = 5000, verbose = TRUE)
})
dags

path <- "./alzheimers44770VC_5000.csv"
geneexpr <- read.csv(path, check.names = FALSE)
dat <- sparsebnData(geneexpr, type = "c")
lambdas <- generate.lambdas(sqrt(nrow(geneexpr)), scale = "linear", lambdas.length = 20)
system.time({
  dags <- estimate.dag(dat, lambdas = lambdas, edge.threshold = 5000, verbose = TRUE)
})
dags

path <- "./alzheimers44771CR_5000.csv"
geneexpr <- read.csv(path, check.names = FALSE)
dat <- sparsebnData(geneexpr, type = "c")
lambdas <- generate.lambdas(sqrt(nrow(geneexpr)), scale = "linear", lambdas.length = 20)
system.time({
  dags <- estimate.dag(dat, lambdas = lambdas, edge.threshold = 5000, verbose = TRUE)
})
dags

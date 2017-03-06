#---------------SECTION 4.1: SPEED AND SCALABILITY IMPROVEMENTS----------------#
#
# This script reproduces the results in Section 4.1 (Figure 2)
#

library("sparsebn")
library("bnlearn")
library("Matrix")
library("mvtnorm")

set.seed(1) # for reproducibility

gen_block_dag <- function(nblocks = 2){
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

gen_gaussian_data <- function(coefs, nn = 50){
    nnode <- ncol(coefs)

    ### Randomly flip half the signs
    nedge <- sum(coefs>0)
    prop.flip <- 0.5
    flipped <- sample(which(coefs>0), floor(prop.flip * nedge))
    coefs[flipped] <- -1 * coefs[flipped]

    ### Generate random data
    vars <- Matrix::Diagonal(n = nnode, x = rep(1, nnode))
    id <- vars
    covMat <- t(solve(id - coefs)) %*% vars %*% solve(id - coefs)
    gaussian.data <- rmvnorm(n = nn, sigma = as.matrix(covMat))

    sparsebnData(gaussian.data, type = "c")
}

run_scalability_test <- function(nblocks, nn = 50, algs = c("ccdr", "pc", "mmhc")){
    stopifnot(is.list(nblocks))

    timing <- vector("list", length = length(nblocks))
    for(k in seq_along(nblocks)){
        dat <- gen_gaussian_data(coefs = gen_block_dag(nblocks[[k]]),
                                 nn = nn)

        times <- rep(NA, 8)
        edges <- rep(NA, 4)
        names(edges) <- c("ccdr", "pc", "mmhc", "hc")

        ### Comparison times
        if("pc" %in% algs){
            library(pcalg)
            times[3] <- proc.time()[3]
            out <- pc(suffStat = list(C = cor(dat$data), n = nrow(dat$data)),
                      indepTest = gaussCItest, ## indep.test: partial correlations
                      alpha=0.01, labels = names(dat$data), verbose = FALSE)
            times[4] <- proc.time()[3]

            edges["pc"] <- graph::numEdges(out@graph)
        }

        if("mmhc" %in% algs){
            library(bnlearn)
            times[5] <- proc.time()[3]
            out <- mmhc(dat$data)
            times[6] <- proc.time()[3]

            edges["mmhc"] <- nrow(out$arcs)
        }

        if("hc" %in% algs){
            times[7] <- proc.time()[3]
            out <- hc(dat$data)
            times[8] <- proc.time()[3]

            edges["hc"] <- nrow(out$arcs)
        }

        # Set edge threshold for CCDr
        threshold <- NA
        threshold <- max(edges, na.rm = T)
        if(is.na(threshold)) threshold <- 5000

        if("ccdr" %in% algs){
            library(sparsebn)
            times[1] <- proc.time()[3]
            dags.out <- estimate.dag(data = dat,
                                     lambdas.length = 20,
                                     edge.threshold = threshold,
                                     verbose = FALSE)
            times[2] <- proc.time()[3]

            edges["ccdr"] <- max(num.edges(dags.out))
        }

        timing[[k]] <- list(nblocks = nblocks[[k]],
                            pp = ncol(dat$data),
                            nn = nrow(dat$data),
                            ccdr = times[2] - times[1],
                            pc = times[4] - times[3],
                            mmhc = times[6] - times[5],
                            hc = times[8] - times[7],
                            threshold = threshold,
                            ccdr = edges["ccdr"],
                            pc = edges["pc"],
                            mmhc = edges["mmhc"],
                            hc = edges["hc"])
        print(do.call("rbind", timing))
    }

    do.call("rbind", timing)
}

### Run simulation (takes ~10 minutes)
pathfinder <- run_scalability_test(as.list(1:10), 50)

### Reproduce Figure 2
dat <- pathfinder[, c(2, 4, 5, 6)]
matplot(dat[,1], dat[,-1],
        pch = c("C", "P", "M"),
        xlab = "p (number of nodes)",
        ylab = "Time (s)",
        type = "b",
        col = c("black", "blue", "green")
        )

### Section 1 --------------------------------------------------------
library(sparsebn)
data(pathfinder)
data <- sparsebnData(pathfinder[["data"]], type = "continuous")
dags <- estimate.dag(data)
plotDAG(dags)

# Re-run with timer
system.time({
    dags <- estimate.dag(data)
})

# Figure 1: Plot first 4 non-trivial DAGs in solution path
# Note: First estimate is always the null graph, and is omitted from plots by default
plotDAG(dags[1:5])

### Section 5 --------------------------------------------------------
library(sparsebn)
data(cytometryContinuous)
names(cytometryContinuous)

### Section 5.1 --------------------------------------------------------
cyto.data <- sparsebnData(cytometryContinuous[["data"]],
                          type = "continuous",
                          ivn = cytometryContinuous[["ivn"]])
cyto.data

### Section 5.2 --------------------------------------------------------
cyto.learn <- estimate.dag(cyto.data)
cyto.learn

cyto.learn[[1]]
cyto.learn[[3]]
print(cyto.learn[[2]], maxsize = 15)
get.adjacency.matrix(cyto.learn[[3]])
show.parents(cyto.learn[[3]], c("raf", "pip2"))

cyto.learn50 <- estimate.dag(data = cyto.data,
                             lambdas.length = 50)

### Use a linear scale
cyto.lambdas <- generate.lambdas(lambda.max = 10,
                                 lambdas.ratio = 0.001,
                                 lambdas.length = 10,
                                 scale = "linear")
cyto.lambdas

### Use a log scale
cyto.lambdas <- generate.lambdas(lambda.max = 10,
                                 lambdas.ratio = 0.001,
                                 lambdas.length = 10,
                                 scale = "log")
cyto.lambdas

cyto.learn.log <- estimate.dag(data = cyto.data,
                               lambdas = cyto.lambdas)

### Section 5.3 --------------------------------------------------------
cyto.param <- estimate.parameters(cyto.learn, data = cyto.data)

cyto.param[[3]]

### Section 5.4 --------------------------------------------------------
select(cyto.learn, edges = 8)  # exact match returned
select(cyto.learn, edges = 10) # closest match returned

select(cyto.learn, lambda = 41.75) # closest match returned
select(cyto.learn, lambda = 41.7)  # same output as previous line

select(cyto.learn, index = 4) # exact match returned
cyto.learn[[4]]               # same output as previous line

select.parameter(cyto.learn, cyto.data)

### Section 5.5 --------------------------------------------------------
getPlotPackage()

setPlotPackage("network")
getPlotPackage()

setPlotPackage("graph")
getPlotPackage()

setPlotPackage("igraph")
getPlotPackage()

plot(cyto.learn)

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
plot(cyto.learn[[7]],
     layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
                              igraph::in_circle()),
     vertex.label = get.nodes(cyto.learn),
     vertex.size = 30,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.5
)

# Figure 6
par(mfrow=c(1,3))
setPlotPackage("igraph")
plot(cyto.learn[[7]])
setPlotPackage("network")
plot(cyto.learn[[7]])
setPlotPackage("graph")
plot(cyto.learn[[7]])

setPlotPackage("igraph") # reset to default (not in paper)

### Section 5.6 --------------------------------------------------------
names(cyto.data)
names(cyto.learn)
names(cyto.learn[[1]])

### Section 6.1 --------------------------------------------------------
library(sparsebn)
data(cytometryDiscrete)

cyto.data <- sparsebnData(cytometryDiscrete[["data"]],
                                 type = "discrete",
                                 ivn = cytometryDiscrete[["ivn"]])
cyto.data

cyto.learn <- estimate.dag(data = cyto.data)
cyto.learn

plot(select(cyto.learn, edges = 15),
     layout = igraph::layout_(to_igraph(cytometryContinuous[["dag"]]),
                              igraph::in_circle()),
     vertex.label = get.nodes(cyto.learn),
     vertex.size = 30,
     vertex.label.color = gray(0),
     vertex.color = gray(0.9),
     edge.color = gray(0),
     edge.arrow.size = 0.5
)

cyto.param <- estimate.parameters(cyto.learn, data = cyto.data)
cyto.param[[5]][["raf"]]

### Section 6.2 --------------------------------------------------------
data(pathfinder)
dat <- sparsebnData(pathfinder[["data"]], type = "c")
nn <- num.samples(dat)
lambdas <- generate.lambdas(sqrt(nn), 0.05,
                            lambdas.length = 50,
                            scale = "linear")
dags <- estimate.dag(data = dat,
                     lambdas = lambdas,
                     edge.threshold = 500,
                     verbose = FALSE)
dags

dat <- sparsebnData(pathfinder[["data"]][1:50, ], type = "c")
nn <- num.samples(dat)
lambdas <- generate.lambdas(sqrt(nn), 0.05,
                            lambdas.length = 50,
                            scale = "linear")
dags <- estimate.dag(data = dat,
                     lambdas = lambdas,
                     edge.threshold = 500,
                     verbose = FALSE)
dags

### Section 6.3 --------------------------------------------------------
set.seed(1) # for reproducibility
nnode <- nedge <- 5000
coefs <- random.dag(nnode, nedge)
id <- Matrix::Diagonal(n = nnode) # identity matrix
vars <- id
covMat <- t(solve(id - coefs)) %*% vars %*% solve(id - coefs)

nobs <- 50
gaussian.data <- mvtnorm::rmvnorm(n = nobs, sigma = as.matrix(covMat))
dat <- sparsebnData(gaussian.data, type = "c")

lambdas <- generate.lambdas(sqrt(nobs), scale = "linear", lambdas.length = 20)
system.time({
    dags <- estimate.dag(dat, lambdas = lambdas, edge.threshold = 5000)
})
dags

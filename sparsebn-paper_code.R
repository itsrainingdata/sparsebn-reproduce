### Section 1 --------------------------------------------------------
library(sparsebn)
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

### Section 4.6 --------------------------------------------------------
install.packages("sparsebn")
devtools::install_github(c("itsrainingdata/sparsebnUtils/dev",
                           "itsrainingdata/ccdrAlgorithm/dev",
                           "gujyjean/discretecdAlgorithm/dev",
                           "itsrainingdata/sparsebn/dev"))

### Section 5 --------------------------------------------------------
library("sparsebn")
data("cytometryContinuous")
names("cytometryContinuous")

### Section 5.1 --------------------------------------------------------
cyto.raw <- cytometryContinuous$data
cyto.ivn <- cytometryContinuous$ivn
cyto.data <- sparsebnData(cyto.raw,
                          type = "continuous",
                          ivn = cyto.ivn)
print(cyto.data)
summary(cyto.data)

### Section 5.2 --------------------------------------------------------
cyto.learn <- estimate.dag(cyto.data)
print(cyto.learn)
summary(cyto.learn)

estimate.dag(cyto.data,
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

estimate.dag(cyto.data,
             lambdas = cyto.lambdas)

### Section 5.3 --------------------------------------------------------
whitelist <- matrix(c("pip3", "pip2"), nrow = 1)
estimate.dag(cyto.data,
             whitelist = whitelist)

whitelist

blacklist <- rbind(c("raf", "jnk"),
                   c("jnk", "raf"))
estimate.dag(cyto.data,
             blacklist = blacklist)

blacklist <- specify.prior(roots = "pip3", leaves = "jnk", nodes = names(cyto.data$data))
estimate.dag(cyto.data,
             blacklist = blacklist)

### Section 5.4 --------------------------------------------------------
print(cyto.learn[[1]])
sumary(cyto.learn[[1]])

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
select(cyto.learn, edges = 10) # closest match returned

select(cyto.learn, lambda = 41.75) # closest match returned
select(cyto.learn, lambda = 41.7)  # same output as previous line

select(cyto.learn, index = 4) # exact match returned
cyto.learn[[4]]               # same output as previous line

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
plot(cyto.learn[[selected.lambda]],
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
par(mfrow=c(1,3))
setPlotPackage("igraph")
plot(cyto.learn[[selected.lambda]])
setPlotPackage("network")
plot(cyto.learn[[selected.lambda]])
setPlotPackage("graph")
plot(cyto.learn[[selected.lambda]])

setPlotPackage("igraph") # reset to default

### The code below requires the RCy3 package to be installed
### along with the Cytoscape application (http://www.cytoscape.org/).
###
### Cytoscape must be installed separately.
source("https://bioconductor.org/biocLite.R")
biocLite("RCy3")

openCytoscape(cytometryContinuous$dag)
openCytoscape(cyto.learn[[selected.lambda]])

### Section 6.1 --------------------------------------------------------
library("sparsebn")
data("cytometryDiscrete")

cyto.data <- sparsebnData(cytometryDiscrete$data,
                          type = "discrete",
                          ivn = cytometryDiscrete$ivn)
cyto.data

cyto.learn <- estimate.dag(cyto.data)
cyto.learn

plot(select(cyto.learn, edges = 19),
     layout = igraph::layout_(to_igraph(cytometryContinuous$dag),
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
dat <- sparsebnData(pathfinder$data, type = "c")

nn <- num.samples(dat)
lambdas <- generate.lambdas(sqrt(nn), 0.05,
                            lambdas.length = 50,
                            scale = "linear")
dags <- estimate.dag(data = dat,
                     lambdas = lambdas,
                     edge.threshold = 500,
                     verbose = FALSE)
dags

dat <- sparsebnData(pathfinder$data[1:50, ], type = "c")
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

#
# See sparsebn-alzheimers.R online: https://github.com/itsrainingdata/sparsebn-reproduce/blob/master/sparsebn-alzheimers.R
# The raw data can be downloaded from this repo as well (required for script)
#

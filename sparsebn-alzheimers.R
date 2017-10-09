#
# Reproduction code for Section 6.3 of the paper "Learning Large-Scale
#  Bayesian Networks with the sparsebn Package"
#

### This code snippet can be run on any of the three datasets below
path <- "./alzheimers44768PFC_5000.csv"
# path <- "./alzheimers44770VC_5000.csv"
# path <- "./alzheimers44771CR_5000.csv"
geneexpr <- read.csv(path, check.names = FALSE)

library(sparsebn)
dat <- sparsebnData(geneexpr, type = "c")
lambdas <- generate.lambdas(sqrt(nrow(geneexpr)), scale = "linear", lambdas.length = 20)
system.time({
    dags <- estimate.dag(dat,
                         lambdas = lambdas,
                         edge.threshold = 5000,
                         verbose = TRUE)
})

dags

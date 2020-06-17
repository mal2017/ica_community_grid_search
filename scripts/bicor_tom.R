library(WGCNA)
library(feather)
library(tidyverse)

threads <- snakemake@threads[[1]]

X <- read_feather(snakemake@input[[1]])

X <- column_to_rownames(X,"index")

# testing
#X <- X[,1:1000]
#threads <- 4

# default - the intermediate correlation matrix/adjacency is assumed unsigned
X.tom.abs <- WGCNA::TOMsimilarityFromExpr(X,
  corType="bicor",
  pearsonFallback='individual',
  verbose=3,
  nThreads = threads)

# signed - the intermediate correlation matrix/adjacency is assumed signed
X.tom.signed <- WGCNA::TOMsimilarityFromExpr(X,
    corType="bicor",
    pearsonFallback='individual',
    networkType="signed",
    verbose=3,
    nThreads = threads)

colnames(X.tom.abs) <- colnames(X)
rownames(X.tom.abs) <- colnames(X)

colnames(X.tom.signed) <- colnames(X)
rownames(X.tom.signed) <- colnames(X)

res.abs <- as_tibble(X.tom.abs, rownames = "index")
res.signed <- as_tibble(X.tom.signed, rownames = "index")

write_feather(res.abs, snakemake@output[["abs"]])
write_feather(res.signed, snakemake@output[["signed"]])

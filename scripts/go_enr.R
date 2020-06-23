library(tidyverse)
library(topGO)
library(org.Dm.eg.db)

# qval_df <- read_csv("test/ica_20comps_rep1_qvalues.csv.gz")
qval_df <- read_csv(snakemake@input[[1]])
ofl <- snakemake@output[[1]]
qval <- snakemake@params[["qval"]]
ont <- snakemake@wildcards[["ont"]]


# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

ag <- qval_df %>% dplyr::select(module,X1,qval) %>%
  mutate(module=as.numeric(module)) %>%
  arrange(module) %>%
  split(.,.$module) %>%
  map(dplyr::select,X1,qval) %>%
  map(deframe)


# func to select significant
#selection <- function(allScore){ return(allScore < snakemake@params[["qval"]])}
selection <- function(allScore){ return(allScore < qval)}

# ontology gene mapping
allGO2genes <- annFUN.org(whichOnto=ont, mapping="org.Dm.eg.db", ID="ensembl")


res <- ag %>%
  map_df(.f=function(x) {
  GOdata <- new("topGOdata",
                ontology=ont,
                allGenes=x,
                annot=annFUN.GO2genes,
                GO2genes=allGO2genes,
                geneSel=selection,
                nodeSize=10)

  result <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

  allRes <- GenTable(GOdata, weight01 = result, topNodes = 10)

  allRes %>% as_tibble() %>%
    mutate(weight01= as.numeric(weight01)) %>%
    arrange(weight01)

},.id = "cluster")

res %>%
  mutate(score = -log10(weight01)) %>%
  write_csv(ofl)



# cl_df <- json_cl %>% enframe(name = "cluster",value = 'genes')
#
# res <- cl_df %>%
#   head(3) %>%
#   mutate(enr=map(genes,.f=function(x) {
#     ego <- enrichGO(gene          = unlist(x),
#                     keyType="FLYBASE",
#                     ont="CC",
#                     OrgDb         = org.Dm.eg.db,
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.01,
#                     qvalueCutoff  = 0.05,
#                     readable      = TRUE)
#     ego <- gofilter(ego,level=7)
#     ego2 <- simplify(ego)
#     as_tibble(ego2)
#   }))
#
#
# res %>% dplyr::select(cluster,enr) %>%
#   deframe() %>%
#   bind_rows(.id="cluster") %>%
#   group_by(cluster) %>% DT::datatable()
# top_n(n = 10, wt = -log10(p.adjust)) %>% DT::datatable()

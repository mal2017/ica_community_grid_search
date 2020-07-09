library(tidyverse)
#library(jsonlite)
library(topGO)

allGO2genesMF <- annFUN.org(whichOnto="MF", mapping="org.Dm.eg.db", ID="ensembl")
allGO2genesCC <- annFUN.org(whichOnto="CC", mapping="org.Dm.eg.db", ID="ensembl")
allGO2genesBP <- annFUN.org(whichOnto="BP", mapping="org.Dm.eg.db", ID="ensembl")

goterms <- list(MF=allGO2genesMF,CC=allGO2genesCC,BP=allGO2genesBP) %>%
  map(names)

#fls <- Sys.glob("~/amarel/200701_ica_grid_full/metrics/enr/enr_fdr*_*comps_rep*_0.0*qval_*.csv")

fls <- snakemake@input

df <- fls %>%
  tibble(fl = .) %>%
  mutate(rep = str_extract(fl,"(?<=rep).+(?=_0.0\\d{1,2}qval)")) %>%
  mutate(ica = paste0("fdr",str_extract(fl,"(?<=fdr)\\d{1}"))) %>%
  mutate(qval = str_extract(fl,"(?<=rep\\d{1}_).+(?=qval)")) %>%
  mutate(comps = str_extract(fl,"(?<=enr_fdr\\d{1,2}_)\\d+(?=comps)")) %>%
  mutate(ont = str_extract(fl,"(?<=qval_).+(?=\\.csv)")) %>%
  #sample_n(100) %>%
  mutate(df = map(fl,.f=function(x) read_csv(x, col_types="cccddddd") %>% group_by(cluster) %>% top_n(5,score))) %>%
  unnest(df) %>%
  dplyr::select(-fl) %>%
  mutate_at(vars(c('comps','qval')),as.numeric)


##
# Percentage of GO terms discovered
##
go_coverage_df <- df %>% group_by(rep,qval,comps,ont,ica) %>%
  summarize(cov=sum(goterms[[head(ont,1)]] %in% GO.ID)/length(goterms[[head(ont,1)]])) %>%
  ungroup() %>%
  group_by(qval,comps,ont,ica) %>%
  summarize(cov=mean(cov)) %>%
  ungroup()

##
# significant terms per module
##

signif_terms_per_mod <- df %>%
  filter(weight01 < 0.05) %>%
  group_by(rep,qval,comps,ont,ica, cluster) %>%
  summarize(n=n()) %>%
  group_by(rep,ont, ica, qval,comps) %>%
  summarize(n=mean(n)) %>%
  group_by(ont, ica, qval,comps) %>%
  summarize(signif_terms_per_mod=mean(n)) %>% arrange(-signif_terms_per_mod) %>%
  ungroup()


##
# percentage mods with at least 1 signif enrichment
##

pct_enr_mods <- df %>%
  group_by(rep,qval,comps,ont,ica, cluster) %>%
  summarize(n=sum(weight01 < 0.05)/n()) %>% 
  group_by(rep,ont, ica, qval,comps) %>%
  summarize(n=mean(n)) %>%
  group_by(ont, ica, qval,comps) %>%
  summarize(pct_enr_mods=mean(n)) %>% arrange(-pct_enr_mods) %>%
  ungroup() 

##
# for each trial for each module, percentage of terms that are uniquely assigned to that module
# THis is then averaged across modules and then across replicates
##

mean_unique_terms <- df %>%
  group_by(rep,ica,qval,comps, GO.ID) %>%
  mutate(is_unique_term = n() == 1 & weight01 < 0.05) %>% 
  group_by(rep,ica,qval,comps, ont,cluster) %>%
  top_n(5,score) %>%
  summarise(pct_unique_term = sum(is_unique_term)/n()) %>%
  group_by(ica,qval,comps, ont, rep) %>%
  summarise(pct_unique_term = mean(pct_unique_term)) %>%
  group_by(ica,qval,comps, ont) %>%
  summarise(pct_unique_term_per_mod = mean(pct_unique_term)) %>%
  arrange(desc(pct_unique_term_per_mod))


##
# Bring all these metrics togehter
##

df2 <- full_join(go_coverage_df,signif_terms_per_mod) %>%
  full_join(mean_unique_terms) %>%
  full_join(pct_enr_mods)

write_csv(df2,snakemake@output[[1]])
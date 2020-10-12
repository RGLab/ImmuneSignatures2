# Pre-vaccination Manuscript data generation
library(tidyverse)

gmxFile <- "data-raw/msigdb.v7.1.symbols.gmt"
# read hallmark
cNames <- count_fields(file = gmxFile, tokenizer = tokenizer_tsv()) %>%
  max() %>%
  seq(from = 1) %>%
  paste0("X", .)
hallmark_geneset <- read_tsv(file = gmxFile, col_names = cNames) %>%
  filter(grepl(pattern = "HALLMARK_", X1)) %>%
  select(-X2) %>%
  pivot_longer(cols = -X1) %>%
  filter(!is.na(value)) %>%
  select(value, X1) %>%
  unstack()
usethis::use_data(hallmark_geneset)

# load BTM (from qusage)
load(file = "data-raw/GeneSets.rda")
btmLS <- stack(BTM.geneSets) %>%
  as.data.frame() %>%
  rowid_to_column()
btmLS$values %>%
  strsplit(split = " /// ") %>%
  setNames(nm = btmLS$rowid) %>%
  stack() %>%
  as.data.frame() %>%
  merge(x = btmLS, by.x = "rowid", by.y = "ind") %>%
  select(values.y, ind) %>%
  unstack() -> btmLS
# append to list of genesets
gmx <- c(gmx, btmLS)

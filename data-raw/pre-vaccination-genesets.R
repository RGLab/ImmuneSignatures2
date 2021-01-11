# Pre-vaccination Manuscript data generation
library(tidyverse)
library(UpdateAnno)
library(data.table)
library(here)

gmxFile <- here("data-raw", "msigdb.v7.1.symbols.gmt")
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

# load BTM (from qusage)
load(file = here("data-raw", "GeneSets.rda"))
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
all_genesets <- c(hallmark_geneset, btmLS)

# Update to latest anno
if (packageVersion("UpdateAnno") != "3.9.0") stop("Install UpdateAnno v3.9.0!")
data("hgncAlias2Symbol")
# Copied from updateDataWithLatestHGNCMap vignette from UpdateAnno
updateModuleSymbols <- function(module, hgncAlias2Symbol){

  # read in as DT
  modList <- get(module)
  modDT <- rbindlist(lapply(modList, data.table), idcol = "module")
  setnames(modDT, "V1", "ALIAS")

  # update gene symbols
  modDT[hgncAlias2Symbol, SYMBOL := SYMBOL, on = c(ALIAS = "ALIAS")]
  setnames(modDT, "module", "pathway")
  modDT[, ALIAS := NULL]
  modDT <- modDT[ !is.na(modDT$SYMBOL) ]
  modDT <- unique(modDT) # MUST DE-DUPE so that gmt file can be created, checked by GeneSet()

  # save with correct name
  upMod <- plyr::dlply(modDT, 1, function(x){ x$SYMBOL }) # convert to list of lists

  return(upMod)
}
all_genesets <- updateModuleSymbols("all_genesets", hgncAlias2Symbol)
btm_list <- updateModuleSymbols("btmLS", hgncAlias2Symbol)
usethis::use_data(all_genesets, overwrite = TRUE)
usethis::use_data(btm_list, overwrite = TRUE)

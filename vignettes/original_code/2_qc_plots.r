## QC PLOTS 
suppressMessages(library(tidyverse)) 
suppressMessages(library(arrayQualityMetrics)) 
suppressMessages(library(here)) 
options(stringsAsFactors = FALSE)

# setup datapath 
datapath = here("QC/generated_data/")
figpath = here("QC/figures/")
# dir.create(datapath) ; dir.create(figpath)

# make machine readable data frame of results 
qc_list = readRDS("QC/generated_data/qclist.rds")
qc.table.names = c( "abs_mean_distance","KS_signal_intensity","D_MAplot")
qc_df = bind_rows(qc_list)
colnames(qc_df)[3:5] = qc.table.names
write_delim(qc_df, path = paste0(datapath,"qc_metadata.txt"), delim = "\t") 

# change "x" vlaues to numeric binary 
change_x = function(x){ if_else(x == "x", 1 ,false =  0)}
qc_df = qc_df %>%  mutate_at(.vars = qc.table.names, change_x)

# qc heatmap 
cu = pals::brewer.blues(n = 3); cu[1] = "#FFFFFF"
pheatmap::pheatmap(qc_df[, 3:5], 
                   cluster_cols = T, cluster_rows = T,
                   color = cu,width = 2, show_rownames = FALSE,
                   filename = paste0(figpath,"outlier_heatmap.png"))
pheatmap::pheatmap(qc_df[, 3:5], 
                   cluster_cols = T, cluster_rows = T,
                   color = cu,width = 2, show_rownames = FALSE,
                   filename = paste0(figpath,"outlier_heatmap.pdf"))

# tidy
dfp = qc_df %>%
  select(qc.table.names,  matrix, age_imputed, gender, timepoint = time_post_last_vax, uid) %>% 
  gather(metric, value, qc.table.names[1]:qc.table.names[length(qc.table.names)]) %>% 
  group_by( matrix, age_imputed, gender, timepoint, uid, value ) %>% 
  summarize(qcsum = sum(value)) %>% 
  ungroup() %>% 
  mutate(age_imputed = as.numeric(age_imputed)) %>% 
  mutate(SDY = str_sub(matrix, 1,7))

p = ggplot(dfp, aes(x = SDY, y = qcsum)) +
  theme_bw() + 
  geom_jitter(shape = 21, size = 0.5, width = 0.3, height = 0.3, fill = "blue3") +
  coord_flip() + 
  geom_hline(yintercept = c(0.5, 1.5, 2.5), linetype = "dashed") +
  # xlab("SDY unique Matrix ") + 
  scale_y_continuous(breaks = c(0,1,2,3)) + 
  # ylim(c(0,3)) + 
  ylab("number of QC \n metrics not passed ")  +
  theme(text = element_text(face = "bold"))
ggsave(p,filename =  paste0(figpath, "array_quality.png"), width = 3, height = 7.5)
ggsave(p,filename =  paste0(figpath, "array_quality.pdf"), width = 3, height = 5)

# venn diagram 
plt = 
  qc_df%>% 
  select(qc.table.names) %>% 
  mutate_if(is.integer, as.numeric) %>% 
  mutate_if(is.logical, as.numeric) 
pdf(file = paste0(figpath, "venn_outlier.pdf"), width = 5, height = 5)
venn::venn(plt, zcolor = "style",
           borders = FALSE, 
           size = 10, 
           cexil = 1.1,
           cexsn = 0.6, 
           ggplot = TRUE) 
dev.off()
#pheatmap::pheatmap(as.matrix(table(dfp$matrix, dfp$qcsum) ),color = cu,
#                   display_numbers = TRUE, number_format = "%.0f", treeheight_row = 0, treeheight_col = 0) 



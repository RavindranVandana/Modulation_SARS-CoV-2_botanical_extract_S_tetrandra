library("limma")
library("dplyr")
####### Code written by Vandana Ravindran University of Oslo ###############
#### function to calculate fold-change using Limma 

get_fold_change <- function(expr_file, annotation_file, out_folder) {
  expr_df <- read.csv(expr_file)
  annot_df <- read.csv(annotation_file)
  
  #### creating the matrix for input
  obj <- expr_df[,-1] %>% as.matrix()
  rownames(obj) <- expr_df[[1]]
  group <- factor(annot_df[, 2])
  
  #### creating design matrix
  design <- model.matrix( ~ 0 + group)
  cont.matrix <- makeContrasts(
    trtvsctrl = (groupcase - groupcontrol),
    levels = design)
  ##### model fitting
  fit = lmFit(log2(obj), design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  fit2e <- eBayes(fit2)
  #### calculating foldchnage 
  fold_change <- topTable(fit2e, number = nrow(obj)) %>%
    mutate(
      logabsFC = abs(logFC),
      FC = 2 ^ (logabsFC),
      rankFC = rank(-FC)
    ) %>% tibble::rownames_to_column(var = "gene_name")
  ##### shortlist ranked genes
  shortlisted <- fold_change %>% 
    filter(logabsFC >= 1) %>% 
    slice_min(rankFC, n = 100) %>%
    mutate(regulation = if_else(logFC <= 0, "down-regulated", "up-regulated"))
  ###### write the file/ creating output 
  fname <- grep("csv", strsplit(expr_file, "/")[[1]], value = TRUE)
  fname <- strsplit(fname, "\\.")[[1]][1]
  write.csv(fold_change, file = file.path(out_folder, paste0(fname, "-foldchange.csv")),row.names = FALSE)
  write.csv(shortlisted, file = file.path(out_folder, paste0(fname, "-shortlisted.csv")),row.names=FALSE)
}

#### mention your path 

path <- "~/Library/folder1/yourpath"

expr_files <- dir(file.path(path, "expression_data"), full.names = TRUE)
annot_files <- dir(file.path(path, "annotation_data"), full.names = TRUE)

for(idx in seq_along(expr_files)) {
  get_fold_change(expr_files[idx], annot_files[idx], "foldchange_output")
}


##### if you want to read each file individually  

# expr_file <- file.path(path, "expression_data", "DMSO_EIDD.csv")
# annot_file <- file.path(path, "annotation_data", "annotation_EIDD.csv")

# get_fold_change(expr_file, annot_file, "foldchange_output")

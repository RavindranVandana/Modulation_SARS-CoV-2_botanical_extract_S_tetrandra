library(tidyverse)

path <- "~/Library/yourpath"
protein_dataset <- read.csv("~/Library/yourpath/protein_dataset.csv")
shortlisted_files <- dir(file.path(path, "foldchange_output"), full.names = TRUE,pattern = "shortlisted")

get_proteinclass <- function(shortlisted, proteindata = protein_dataset) {
    shortlist_df <- read.csv(shortlisted)
    proteinclass <- merge(
      x = shortlist_df,
      y = proteindata,
      by.x = "gene_name",
      by.y = "From",
      all.x = TRUE
    )
    proteinclass %>% 
      distinct(gene_name, l1) %>% 
      count(l1) %>% 
      filter(!is.na(l1) & !grepl("Unclass", l1)) %>% 
      mutate(Group = gsub(".*0.*? (.*)-shortlisted\\.csv", "\\1",  shortlisted))
  }

write.csv(temp,"temp.csv",row.names = FALSE)

protein_class<-purrr::map_df(shortlisted_files, get_proteinclass)

protein_class<-protein_class %>% mutate(Compound = recode(Group, DMSO_EIDD = 'MPV', DMSO_fang = 'Fangchinoline', DMSO_frac =  'Fraction18' ,DMSO_tetra="Tetrandrine"))

ggplot(data = protein_class, aes(y = l1, x = n, fill = l1)) +
  geom_bar(stat = "identity", position = "dodge") +
  
scale_fill_viridis_d() +
  labs(x = "Number of genes",
       y = "Protein class",
       fill = "Compound") +
  theme_minimal(base_size = 15) + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor= element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 15, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 15, colour = "black")
  )+
  scale_x_continuous(breaks = seq(0,35,5),expand = expansion())+facet_wrap(~Compound)

##################### heat map for fold-change ######################################

########## heat map of all viral genes############

########## creating the file #############
library(readxl)
viral_genes <-as.data.frame(read_excel("viral_genes.xlsx")) %>% mutate(Symbol=toupper(trimws(Viral_genes))) 
viral_genes<-viral_genes%>% mutate(Symbol=toupper(trimws(Viral_genes))) 
foldchange_files <- dir(file.path(path, "foldchange_output"), full.names = TRUE,pattern = "foldchange")



get_viralgenes <- function(foldchange, viraldata = viral_genes) {
  foldchange_df <- read.csv(foldchange) %>% 
    mutate(id = toupper(trimws(gene_name)))
  viralgenes <- merge(
    x = foldchange_df,
    y = viraldata,
    by.x = "id",
    by.y = "Symbol",
    all.y = TRUE
  ) %>% dplyr::select(id, logFC) %>% 
    mutate(Group = gsub(".*0.*? (.*)-foldchange\\.csv", "\\1",  foldchange))
  
  return(viralgenes)
}

viralgenes<-purrr::map_df(foldchange_files, get_viralgenes)
viralgenes <-
  viralgenes %>% mutate(
    Compound = recode(
      Group,
      DMSO_EIDD = 'MPV',
      DMSO_fang = 'Fangchinoline',
      DMSO_frac =  'Fraction18' ,
      DMSO_tetra = "Tetrandrine"
    ),
    id=recode(id,
    `SARS-COV-2_E`="E",
    `SARS-COV-2_M`='M',
    `SARS-COV-2_N`='N',
    `SARS-COV-2_S`='S',
    `SARS-COV-2_ORF1AB`='ORF1ab',
    `SARS-COV-2_ORF3A`='ORF3a',
    `SARS-COV-2_ORF7A`='ORF7a',
    `SARS-COV-2_ORF8`='ORF8')
  )


write.csv(viralgenes,"viral_data.csv",row.names = FALSE)      

##### heatmap

viralgenes <- viralgenes %>% 
  mutate(Compound = factor(Compound)) %>% 
  mutate(Compound = reorder(Compound, logFC)) %>% 
  mutate(Compound = forcats::fct_relevel(Compound, "MPV", after = Inf)) %>% 
  mutate(id=forcats::fct_relevel(id, "S", after = 0))

viralgenes %>% 
  ggplot(aes(Compound, id, fill = logFC)) +
  geom_tile(color = "whitesmoke", size = 0.1) +
  coord_equal() +
  scale_fill_viridis_c(option = "D") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text = element_text(size = 15, color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 15,colour = 'black')
  ) + labs(x =NULL, y = NULL, fill = "Log(FC)")


##################################### Clusterprofile analysis#####################################

####### Frac data 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

#### get uniprot_ids

uniprot_ids_ncounter <- read_excel("uniprot_ids_ncounter.xlsx")

get_uniprot<-function(data,protein=uniprot_ids_ncounter){
  shortlist_df <- read.csv(data)
    proteinclass <- merge(
    x = shortlist_df,
    y = protein[,c("From","Entry")],
    by.x = "gene_name",
    by.y = "From",
    all.x = TRUE
  )
    return(proteinclass)
}



frac_genes<-get_uniprot(shortlisted_files[3]) %>% 
  filter(Entry !="NA") 
  
frac_genes_convert<-bitr(frac_genes$Entry, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")


frac_data_cmp <-
  inner_join(x = frac_genes_convert,
             y = frac_genes,
             by = c("UNIPROT" = "Entry")) %>% 
  dplyr::select(ENTREZID, logFC, regulation) 

frac_data_cmp$group <- "Fraction18"

######### tetra data 
tetra_genes<-get_uniprot(shortlisted_files[4]) %>% 
  filter(Entry !="NA") 

tetra_gene_convert<-bitr(tetra_genes$Entry, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")

tetra_data_cmp <-
  inner_join(x = tetra_gene_convert,
             y = tetra_genes,
             by = c("UNIPROT" = "Entry")) %>% 
  dplyr::select(ENTREZID, logFC, regulation) 

tetra_data_cmp$group <- "Tetrandrine"


###### fang data 
fang_genes<-get_uniprot(shortlisted_files[2]) %>% filter(Entry !="NA") 
fang_gene_convert<-bitr(fang_genes$Entry, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")

fang_data_cmp <-
  inner_join(x = fang_gene_convert,
             y = fang_genes,
             by = c("UNIPROT" = "Entry")) %>% 
  dplyr::select(ENTREZID, logFC, regulation) 

fang_data_cmp$group <- "Fangchinoline"

####### EIDD data
EIDD_genes<-get_uniprot(shortlisted_files[1]) %>% filter(Entry !="NA") 
EIDD_gene_convert<-bitr(EIDD_genes$Entry, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")


EIDD_data_cmp <-
  inner_join(x = EIDD_gene_convert,
             y = EIDD_genes,
             by = c("UNIPROT" = "Entry")) %>% 
  dplyr::select(ENTREZID, logFC, regulation) 

EIDD_data_cmp$group <- "MPV"



clusterdata<-rbind(frac_data_cmp,tetra_data_cmp,fang_data_cmp,EIDD_data_cmp)


formula_res <-
  compareCluster(
    ENTREZID ~ regulation + group,
    data = clusterdata,
    fun = "enrichGO",
    OrgDb = 'org.Hs.eg.db',
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )



plot<-dotplot(formula_res, x="regulation") + facet_grid(~group)+labs(x="",fill="p-adjusted")

plotkegg<-dotplot(formula_kegg, x="regulation") + facet_grid(~group)+labs(x="",fill="p-adjusted")

plot + scale_color_continuous(name="p-adjusted",type = "viridis")+scale_size_continuous(guide = 'none')

plot_data<-plot$data





########## Heat map #####
complete(plot_data, regulation, group, Description) %>% mutate(group = forcats::fct_relevel(group, "MPV", after = Inf)) %>% 
  ggplot(aes(group, reorder(Description,p.adjust), fill = p.adjust)) +
  geom_tile(color = "whitesmoke", size = 0.1) +
  coord_equal() +
  scale_fill_viridis_c(
    breaks = scales::breaks_log(n = 15),
    option = "D",
    na.value = 'grey',
    direction =-1, 
    trans = "log", 
    labels = scales::label_log(digits = 1),
    guide = guide_colorbar(
      barheight = unit(0.5, "npc")
    )
  ) +
  scale_alpha_continuous(range = c(0.5, 1)) +
  facet_grid(~regulation)+
  labs(x = NULL, #"Regulation", 
       y = NULL, # "GO Biological process",
       fill="p-adjusted") +theme_bw()+
  theme(axis.text.x = element_text(
    size = 15,
    angle = 45,
    vjust = 1,
    hjust = 1
    
  ))  +
  theme(
    axis.ticks = element_blank(),
    legend.text.align = 0,
    legend.position = "right",
    text = element_text(size = 15),
    axis.text.x = element_text(size = 15,colour = "black"),
    axis.text.y = element_text(size = 15,colour = "black"),
    strip.text = element_text(size=15, face="bold", color="black")) +
    scale_y_discrete(
    labels = function(x)
      stringr::str_wrap(stringr::str_to_sentence(x), width = 35),
    expand = expansion()
  ) + scale_x_discrete(expand = expansion())
  






library(DESeq2)
library(tidyverse)
library(readxl)
library(edgeR)
library(ggrepel)
library(fgsea)
library(gridExtra)
library(ComplexHeatmap)
library(testit)
library(stringi)
library(assertable)
library(dplyr)
library(circlize)
library(biomaRt)
library(dplyr)

volcanoPlotWithList <- function(result, fdr, l2fc_cutoff=1, title, gene_list){
  #'
  #' @param result DESeq2 results table with log2Foldchange, padj and Gene_name as columns
  #' @param fdr determines fdr cutoff
  #' @param l2fc_cutoff determines log2Foldchange cutoff
  #' @param title Determines title of the plot
  #'
  
  #'
  assert_colnames(result, colnames=c("padj", "log2FoldChange", "Gene_name"), only_colnames = F, quiet = T)
  ggplot(data = result, mapping = aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(data = result %>% filter (!Gene_name %in% gene_list), color="darkgrey", alpha=0.35)+
    geom_point(data = result %>% filter (Gene_name %in% gene_list), color="navy")+
    geom_text_repel(data = result %>% filter(Gene_name %in% gene_list & -log10(padj) > -log10(fdr)), aes(label=Gene_name), size=2, max.overlaps = Inf)+
    geom_vline(xintercept = -l2fc_cutoff, linetype="dashed")+
    geom_vline(xintercept = l2fc_cutoff, linetype="dashed")+
    geom_hline(yintercept = -log10(fdr), linetype="dashed")+
    theme_bw()+
    ggtitle(title)
}

volcanoPlot.padj.sorted <- function(result, fdr, l2fc_cutoff=1, title){
  #'
  #' @param result DESeq2 results table with log2Foldchange, padj and Gene_name as columns
  #' @param fdr determines fdr cutoff
  #' @param l2fc_cutoff determines log2Foldchange cutoff
  #' @param title Determines title of the plot
  #'
  assert_colnames(result, colnames=c("padj", "log2FoldChange", "Gene_name"), only_colnames = F, quiet = T)
  # Split it in two parts and merge then
  
  # Sort the data frame by padj for upregulated genes
  df_upregulated <- result%>% filter(log2FoldChange > l2fc_cutoff & -log10(padj) > -log10(fdr))
  df_upregulated_sorted <- df_upregulated[order(df_upregulated$padj), ]
  
  # Sort the data frame by padj for downregulated genes
  df_downregulated <- result %>%filter(log2FoldChange < -l2fc_cutoff & -log10(padj) > -log10(fdr))
  df_downregulated_sorted <- df_downregulated[order(df_downregulated$padj), ]
  
  # Subset the top 20 most significant genes for upregulated and downregulated separately
  top_twenty_upregulated <- head(df_upregulated_sorted, 20)
  top_twenty_downregulated <- head(df_downregulated_sorted, 20)
  
  # Create the volcano plot
  volcano_plot <- ggplot(data = result, mapping = aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(data = result %>% filter (abs(log2FoldChange <= l2fc_cutoff) | -log10(padj) <= -log10(fdr)), color="darkgrey")+
    geom_point(data = result %>% filter(log2FoldChange > l2fc_cutoff & -log10(padj) > -log10(fdr)), color="tomato3")+
    geom_point(data = result %>% filter(log2FoldChange < -l2fc_cutoff & -log10(padj) > -log10(fdr)), color="seagreen4")+
    geom_vline(xintercept = -l2fc_cutoff, linetype="dashed")+
    geom_vline(xintercept = l2fc_cutoff, linetype="dashed")+
    geom_hline(yintercept = -log10(fdr), linetype="dashed")+
    theme_bw()+
    ggtitle(title)
  
  # Add labels for the top 20 upregulated genes
  volcano_plot_with_labels <- volcano_plot +
    geom_text_repel(data = top_twenty_upregulated, aes(label = Gene_name), size=1.75, max.overlaps = Inf)
  
  # Add labels for the top 20 downregulated genes
  volcano_plot_with_labels <- volcano_plot_with_labels +
    geom_text_repel(data = top_twenty_downregulated, aes(label = Gene_name), size=1.75,  max.overlaps = Inf)

  # Display the plot
  return(volcano_plot_with_labels)
}



volcanoPlot <- function(result, fdr, l2fc_cutoff=1, title){
  #'
  #' @param result DESeq2 results table with log2Foldchange, padj and Gene_name as columns
  #' @param fdr determines fdr cutoff
  #' @param l2fc_cutoff determines log2Foldchange cutoff
  #' @param title Determines title of the plot
  #'
  
  #'
  assert_colnames(result, colnames=c("padj", "log2FoldChange", "Gene_name"), only_colnames = F, quiet = T)
  ggplot(data = result, mapping = aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(data = result %>% filter (abs(log2FoldChange <= l2fc_cutoff) | -log10(padj) <= -log10(fdr)), color="darkgrey")+
    geom_point(data = result %>% filter(log2FoldChange > l2fc_cutoff & -log10(padj) > -log10(fdr)), color="tomato3")+
    geom_point(data = result %>% filter(log2FoldChange < -l2fc_cutoff & -log10(padj) > -log10(fdr)), color="seagreen4")+
    geom_text_repel(data = result %>% filter(abs(log2FoldChange) > l2fc_cutoff & -log10(padj) > -log10(fdr)), aes(label=Gene_name), size=2)+
    geom_vline(xintercept = -l2fc_cutoff, linetype="dashed")+
    geom_vline(xintercept = l2fc_cutoff, linetype="dashed")+
    geom_hline(yintercept = -log10(fdr), linetype="dashed")+
    theme_bw()+
    ggtitle(title)
}

# Filtered_named_TPM_raw_counts
# Read Aldara Recall data
Aldara_raw_count <- readRDS("./Aldara_counts.rds")
Aldara_raw_count$ensembl_id <- gsub("(ENSG[0-9]+).*","\\1", as.character(Aldara_raw_count$ensembl_id))
dim(Aldara_raw_count) #7 samples; 61533 genes
Aldara_raw_count <- Aldara_raw_count[(rowSums(Aldara_raw_count[,-1])>0),]
any(duplicated(Aldara_raw_count)) #no duplicated rows left

# Metadata Aldara
sample_combi <- read_excel("./sample_matrix.xlsx") 
sample_combi <- sample_combi[order(sample_combi$sample),]
# All samples from Aldara Recall
aldara_samples <- sample_combi %>% 
  dplyr::mutate(dataset = "Aldara")%>% 
  dplyr::rename(sampleID = sample)


# Read in BRAIN
BRAIN_SummarizedExperiment <- readRDS("./dds_highQual_Sexandbatchcorrected_v04.rds")
gene_info <- rowData(BRAIN_SummarizedExperiment)
BRAIN_raw_counts <- counts(BRAIN_SummarizedExperiment, normalized=F)
BRAIN_metadata <- data.frame(colData(BRAIN_SummarizedExperiment))
non_lesional_samples <- BRAIN_metadata %>% 
  dplyr::filter(diag == "non-lesional") %>% 
  dplyr::mutate(dataset = "BRAIN") %>%
  dplyr::select(sampleID, dataset)
BRAIN_healthy <- BRAIN_raw_counts[,non_lesional_samples$sampleID]
BRAIN_healthy <- as.data.frame(BRAIN_healthy) %>%
  rownames_to_column() %>%
  dplyr::rename(ensembl_id = rowname)
any(duplicated(BRAIN_healthy)) #no duplicated rows

# Merge them for joint preprocessing
raw_counts_Aldara_BRAIN <- inner_join(Aldara_raw_count, BRAIN_healthy, by = "ensembl_id") 
raw_counts_Aldara_BRAIN_filt <- raw_counts_Aldara_BRAIN[(rowSums(raw_counts_Aldara_BRAIN[,-1])>0),] #all cols except the ensembl_id one
dim(raw_counts_Aldara_BRAIN_filt) #17653 genes remaining   294 samples     

#add biomaRt gene names

ensembl <- useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl",
                      host="https://may2021.archive.ensembl.org", mirror="www") #https://dec2017.archive.ensembl.org

geneLength_BA <-getBM(attributes=c("ensembl_gene_id",
                                   "ensembl_transcript_id",
                                   "external_gene_name",
                                   "transcript_length",
                                   "transcript_appris",
                                   "hgnc_symbol",
                                   "entrezgene_id",
                                   "description",
                                   "gene_biotype"),
                      filters = "ensembl_gene_id",
                      values = raw_counts_Aldara_BRAIN_filt$ensembl_id,
                      mart=ensembl,
                      useCache = FALSE)


#filter for protein coding genes
geneLength_filt_BA <- geneLength_BA[geneLength_BA$gene_biotype == 'protein_coding', ]

#filter for genes with hgnc_symbol and entrezgene_id 
geneLength_filt_BA <- geneLength_filt_BA[!(is.na(geneLength_filt_BA$entrezgene_id) | geneLength_filt_BA$entrezgene_id==""), ]
geneLength_filt_BA <- geneLength_filt_BA[!(is.na(geneLength_filt_BA$hgnc_symbol) | geneLength_filt_BA$hgnc_symbol==""), ]

#only use principal transcripts, which are annotated as such or average if multiple
principal_BA<- geneLength_filt_BA %>%
  dplyr::filter(grepl("principal", transcript_appris)) %>%
  group_by(ensembl_gene_id,external_gene_name) %>%
  dplyr::summarise(Length= round(mean( transcript_length ))) %>%
  dplyr::rename(GeneID=ensembl_gene_id) 
dim(principal_BA) #17153     3

#add effective transcript length to filtered matrix
length_count_BA <- principal_BA %>%
  dplyr::select(GeneID, Length) %>%
  inner_join(raw_counts_Aldara_BRAIN_filt, by = c("GeneID" = "ensembl_id")) %>%
  mutate(Length=ifelse(Length-500+1>0,Length-500+1,1)) #hiseq 4000 250bp paired end

#calculate TPM -> first correct for transcript length and then for library size
tpm_matrix_BA <-length_count_BA %>%
  pivot_longer(cols = starts_with("MUC"),
               names_to = "sample",
               values_to = "count") %>%
  group_by(sample) %>%
  mutate(rpk = count/(Length/1000) ,
         libSize = sum(rpk),
         tpm = rpk /libSize *1e6 )

#add TPMs 
tpm_BA <- tpm_matrix_BA %>% mutate(minTPM= min(tpm[tpm>0]),
                                   log10TPM = log10(tpm+minTPM)) %>%
  dplyr::select(-libSize,-rpk,-minTPM)
dim(tpm_BA)

#filter tpm > 0 
tpm_filt_BA <- tpm_BA %>% ungroup() %>%
  group_by(GeneID) %>%
  dplyr::filter(tpm>0)  #results in 17,153 genes
filtered_TPM_counts_BA <- tpm_filt_BA %>% dplyr::select(GeneID,sample,count) %>%
  distinct() %>% 
  pivot_wider(names_from = sample,
              values_from = count) 
dim(filtered_TPM_counts_BA) #17153   295

names_BA <- principal_BA %>% dplyr::select(GeneID, external_gene_name)
filtered_named_TPM_counts_BRAIN_Aldara <- inner_join(names_BA, filtered_TPM_counts_BA, by = "GeneID")
###### End of preprocessing ------- ############################################


Aldara_raw_count_filt <- Aldara_raw_count %>% filter(ensembl_id %in% filtered_named_TPM_counts_BRAIN_Aldara$GeneID)
dim(Aldara_raw_count_filt) # 17153 x 8

raw_counts_Aldara_BRAIN_filt <- raw_counts_Aldara_BRAIN_filt %>%
  dplyr::filter(ensembl_id %in% filtered_named_TPM_counts_BRAIN_Aldara$GeneID)
dim(raw_counts_Aldara_BRAIN_filt) # 17153 x 295


genes_filtered_named_counts_Aldara_BRAIN <- raw_counts_Aldara_BRAIN_filt %>% 
  dplyr::select(ensembl_id)
BRAIN_Aldara_counts <- raw_counts_Aldara_BRAIN_filt %>% 
  dplyr::select(!ensembl_id) # 17153 x 295

geneDescription <- getBM(attributes=c("ensembl_gene_id",
                                      "description",
                                      "external_gene_name"),
                         filters = "ensembl_gene_id",
                         values = genes_filtered_named_counts_Aldara_BRAIN, #ensembl IDs
                         mart=ensembl,
                         useCache = FALSE)

#end1_FL:##
end1_FL_Aldara <- aldara_samples %>% 
  dplyr::filter(t_area_p == "end1_FL") %>% 
  dplyr::select(sampleID, dataset)

#end3_FL: Enf of Former lesional Patient
end3_FL_Aldara <- aldara_samples %>% 
  dplyr::filter(t_area_p == "end3_FL") %>% 
  dplyr::select(sampleID, dataset)

#### From here plots for the paper -----------------------------------------------------
#compare Patient1 FL end  to all non-lesional samples of BRAIN
end1_FL_BRAIN_Aldara <- rbind(non_lesional_samples, end1_FL_Aldara)
end1_FL_BRAIN_Aldara <- end1_FL_BRAIN_Aldara[order(end1_FL_BRAIN_Aldara$sampleID),]
rownames(end1_FL_BRAIN_Aldara) <- end1_FL_BRAIN_Aldara$sampleID

# check if normalised BRAIN aldara counts do really contain normalised counts of BRAIN and ALDARA
#order raw counts
sort_end1_FL_counts_Aldara_BRAIN <- BRAIN_Aldara_counts[colnames(BRAIN_Aldara_counts) %in% end1_FL_BRAIN_Aldara$sampleID]
sort_end1_FL_counts_Aldara_BRAIN <- sort_end1_FL_counts_Aldara_BRAIN[,order(colnames(sort_end1_FL_counts_Aldara_BRAIN))]
assert(colnames(sort_end1_FL_counts_Aldara_BRAIN) == end1_FL_BRAIN_Aldara$sampleID)

#create DGEList object
d_end1_FL_lesional <- DGEList(counts = sort_end1_FL_counts_Aldara_BRAIN,
                              genes = genes_filtered_named_counts_Aldara_BRAIN,
                              samples = end1_FL_BRAIN_Aldara)

dataset <- end1_FL_BRAIN_Aldara$dataset
design_nonL <- model.matrix(~0+ dataset)

keep <- filterByExpr(d_end1_FL_lesional, design = design_nonL)
d_end1_FL_lesional <- d_end1_FL_lesional[keep,,keep.lib.sizes=FALSE]

d_end1_FL_lesional <- calcNormFactors(d_end1_FL_lesional)
d_end1_FL_lesional <- estimateDisp(d_end1_FL_lesional, design = design_nonL)
fit_NL_eczema <- glmFit(d_end1_FL_lesional, design_nonL)  

end1_FL_Aldara_vs_nonL_BRAIN <- glmLRT(fit_NL_eczema, contrast = c(1,-1))
end1_FL_Aldara_vs_nonL_BRAIN_all <- as.data.frame(topTags(end1_FL_Aldara_vs_nonL_BRAIN, n = 18000, adjust.method = "BH", sort.by = "PValue", p.value = 1)) %>% 
  dplyr::mutate(DGE = case_when(FDR < 0.05 & logFC > 1 ~ "up",
                                FDR < 0.05 & logFC < 1 ~ "down",
                                T ~ "no")) %>% 
  dplyr::inner_join(geneDescription, by = c("ensembl_id" = "ensembl_gene_id"))
summary(decideTests(end1_FL_Aldara_vs_nonL_BRAIN, adjust.method = "BH", p.value = 0.05)) 
# 39 down 372 up 

# #volcano plot eczema non lesional
end1_FL_Aldara_vs_nonL_BRAIN_all <- end1_FL_Aldara_vs_nonL_BRAIN_all %>% rename(padj =FDR, log2FoldChange=logFC, Gene_name=external_gene_name)
volcanoPlot(end1_FL_Aldara_vs_nonL_BRAIN_all, fdr=0.05, l2fc_cutoff = 1, title = "Lesional eczema Aldara (end patient 1 FL) vs non-lesional BRAIN")

#compare Patient3 end FL to all non-lesional samples of BRAIN
end3_FL_BRAIN_Aldara <- rbind(non_lesional_samples, end3_FL_Aldara)
end3_FL_BRAIN_Aldara <- end3_FL_BRAIN_Aldara[order(end3_FL_BRAIN_Aldara$sampleID),]
rownames(end3_FL_BRAIN_Aldara) <- end3_FL_BRAIN_Aldara$sampleID

#order raw counts
sort_end3_FL_counts_Aldara_BRAIN <- BRAIN_Aldara_counts[colnames(BRAIN_Aldara_counts) %in% end3_FL_BRAIN_Aldara$sampleID]
sort_end3_FL_counts_Aldara_BRAIN <- sort_end3_FL_counts_Aldara_BRAIN[,order(colnames(sort_end3_FL_counts_Aldara_BRAIN))]
assert(colnames(sort_end3_FL_counts_Aldara_BRAIN) == end3_FL_BRAIN_Aldara$sampleID)
#create DGEList object
d_end3_FL_lesional <- DGEList(counts = sort_end3_FL_counts_Aldara_BRAIN,
                              genes = genes_filtered_named_counts_Aldara_BRAIN,
                              samples = end3_FL_BRAIN_Aldara)

dataset <- end3_FL_BRAIN_Aldara$dataset
design_nonL <- model.matrix(~0+ dataset)

keep <- filterByExpr(d_end3_FL_lesional, design = design_nonL)
d_end3_FL_lesional <- d_end3_FL_lesional[keep,,keep.lib.sizes=FALSE]

d_end3_FL_lesional <- calcNormFactors(d_end3_FL_lesional)
d_end3_FL_lesional <- estimateDisp(d_end3_FL_lesional, design=design_nonL)
fit_FL_pso <- glmFit(d_end3_FL_lesional, design_nonL)  

end3_FL_Aldara_vs_nonL_BRAIN <- glmLRT(fit_FL_pso, contrast = c(1,-1))
end3_FL_Aldara_vs_nonL_BRAIN_all <- as.data.frame(topTags(end3_FL_Aldara_vs_nonL_BRAIN, n = 18000, adjust.method = "BH", sort.by = "PValue", p.value = 1)) %>% 
  dplyr::mutate(DGE = case_when(FDR < 0.05 & logFC > 1 ~ "up",
                                FDR < 0.05 & logFC < 1 ~ "down",
                                T ~ "no")) %>% 
  dplyr::inner_join(geneDescription, by = c("ensembl_id" = "ensembl_gene_id"))
summary(decideTests(end3_FL_Aldara_vs_nonL_BRAIN, adjust.method = "BH", p.value = 0.05)) 
# 9 down 75 up 

# #volcano plot eczema non lesional
end3_FL_Aldara_vs_nonL_BRAIN_all <- end3_FL_Aldara_vs_nonL_BRAIN_all %>% rename(padj=FDR, log2FoldChange=logFC, Gene_name = external_gene_name)
volcanoPlot(end3_FL_Aldara_vs_nonL_BRAIN_all, fdr=0.05, l2fc_cutoff = 1, title = "Lesional psoriasis Aldara (end patient 3 FL) vs non-lesional BRAIN")

# End of BRAIN Comparisons -----------------------------------------------------
rownames(Aldara_raw_count_filt) <- Aldara_raw_count_filt$ensembl_id
Aldara_raw_count_filt$ensembl_id <- NULL
dim(Aldara_raw_count_filt)

rownames(sample_combi) <- sample_combi$sample
assert(colnames(Aldara_raw_count_filt) == rownames(sample_combi))

########## 1 vs. 1 comparisons within Aldara ----------------------------#######
BCV = 0.4 # use default BCV for human data as dispersion cannot be estimated due to lack of replicates

# Use Aldara Raw Count Filt
assertthat::are_equal(colnames(Aldara_raw_count_filt), sample_combi$sample)
y <- DGEList(counts = Aldara_raw_count_filt, genes=rownames(Aldara_raw_count_filt), samples = sample_combi)
TPA <- sample_combi$t_area_p
design <- model.matrix(~0+TPA)
y <- calcNormFactors(y)
# Important to take the square of BVC for dispersion
fit <- glmFit(y, design = design, dispersion = BCV**2)
my.contrasts <- makeContrasts(P1end_start = TPAend1_FL - TPAstart1_FL, P3end_start = TPAend3_FL - TPAstart3_FL, P3end_P1end = TPAend3_FL - TPAend1_FL ,levels = design)
 # Different comparisons
Pat1_end_start <- glmLRT(fit, contrast = my.contrasts[,"P1end_start"])
Pat1_end_start <- Pat1_end_start <- topTags(Pat1_end_start, n=length(Aldara_raw_count_filt$MUC20488), adjust.method = "BH")$table

Pat3_FL_end_start <- glmLRT(fit, contrast = my.contrasts[, "P3end_start"])
Pat3_FL_end_start <- topTags(Pat3_FL_end_start, n=length(Aldara_raw_count_filt$MUC20488), adjust.method = "BH")$table

Pat3_FL_end_Pat1_FL_end <- glmLRT(fit, contrast = my.contrasts[, "P3end_P1end"])
Pat3_FL_end_Pat1_FL_end <- topTags(Pat3_FL_end_Pat1_FL_end, n=length(Aldara_raw_count_filt$MUC20488),adjust.method = "BH")$table

Pat1_end_start <- data.frame(merge(Pat1_end_start, gene_info[,c("Gene_name", "ensembl_id")], by.x="genes", by.y="ensembl_id")) %>% dplyr::rename(log2FoldChange=logFC, padj=FDR)
Pat3_FL_end_start <- data.frame(merge(Pat3_FL_end_start, gene_info[,c("Gene_name", "ensembl_id")], by.x="genes", by.y="ensembl_id")) %>% dplyr::rename(log2FoldChange=logFC, padj=FDR)
Pat3_FL_end_Pat1_FL_end <- data.frame(merge(Pat3_FL_end_Pat1_FL_end, gene_info[,c("Gene_name", "ensembl_id")], by.x="genes", by.y="ensembl_id"))%>% dplyr::rename(log2FoldChange=logFC, padj=FDR)

# Figure 2A
volcanoPlot.padj.sorted(Pat1_end_start, fdr = 0.05, title="Patient 1")
volcanoPlot.padj.sorted(Pat3_FL_end_start, fdr = 0.05, title="Patient 3")



# --------------------- ENRICHMENTS HERE --------------------------------------
# LOAD GENE SETS FOR PSORIASIS, ICD AND TRMS
pso_genes <- read_xlsx("./SciTranslMed_2014_Quaranta_Pso genes.xlsx")
pso_genes <- pso_genes %>% dplyr::select(gene)
pso_genes <- pso_genes[!is.na(pso_genes)]
pso_genes <- pso_genes[pso_genes %in% filtered_named_TPM_counts_BRAIN_Aldara$external_gene_name] # 81
ICD_genes <- read_xlsx("./Aldara_gene lists.xlsx", sheet="PL top genes") # 1072
# filter for protein coding genes
ICD_genes_merged <- ICD_genes %>% filter(hgnc_symbol %in% filtered_named_TPM_counts_BRAIN_Aldara$external_gene_name) # 956
# Create a Waterfall plot
unique_ICD_Genes_all <- ICD_genes_merged[!duplicated(ICD_genes_merged$hgnc_symbol),] # 850
unique_ICD_Genes_all$signedP_val <- sign(unique_ICD_Genes_all$fc) * -log10(unique_ICD_Genes_all$pval)
unique_ICD_Genes_all <- unique_ICD_Genes_all[order(unique_ICD_Genes_all$signedP_val, decreasing = T),]
# Repeat pathway analysis with shorter ICD list
shorter_pathways <- list()
shorter_pathways$PSO <- pso_genes
shorter_ICD_names <- unique_ICD_Genes_all%>% filter(abs(signedP_val)>12)
# shorter_ICD_names <- shorter_ICD_names$Gene_name
shorter_ICD_names <- shorter_ICD_names$hgnc_symbol
shorter_pathways$ICD <- shorter_ICD_names

# End Patient1 FL Enrichment
# Figure 2 B HERE------------------------------------------------
# End Patient 1 enrichment
signed_p_values_P1 <- -log10(end1_FL_Aldara_vs_nonL_BRAIN_all$padj) * sign(end1_FL_Aldara_vs_nonL_BRAIN_all$log2FoldChange)
names(signed_p_values_P1) <- end1_FL_Aldara_vs_nonL_BRAIN_all$Gene_name
signed_p_values_P1 <- signed_p_values_P1[order(signed_p_values_P1, decreasing = T)]
signed_p_values_P1 <- signed_p_values_P1[!is.na(signed_p_values_P1)]

# End Patient 3FL Enrichment
# rank the genes according to signedPvalue
signed_P_values_P3 <- -log10(end3_FL_Aldara_vs_nonL_BRAIN_all$padj)*sign(end3_FL_Aldara_vs_nonL_BRAIN_all$log2FoldChange)
names(signed_P_values_P3) <- end3_FL_Aldara_vs_nonL_BRAIN_all$Gene_name
signed_P_values_P3 <- signed_P_values_P3[order(signed_P_values_P3, decreasing = T)]
signed_P_values_P3 <- signed_P_values_P3[!is.na(signed_P_values_P3)]

# Patient1 endFL enrichment Plot
fgsea_shorter1 <- fgsea(shorter_pathways, signed_p_values_P1)
plotEnrichment(shorter_pathways$PSO, signed_p_values_P1)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.5,1.1,0.1), limits = c(-0.5, 1.0))+theme_bw()
plotEnrichment(shorter_pathways$ICD, signed_p_values_P1)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.5,1.1,0.1), limits = c(-0.5, 1.0))+theme_bw()
fgsea_shorter1

# Patient 3 end FL enrichment plot
fgsea_shorter3 <- fgsea(shorter_pathways, signed_P_values_P3)
plotEnrichment(shorter_pathways$PSO, signed_P_values_P3)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.1,1.1,0.1), limits = c(-0.1, 1.0))+theme_bw()
plotEnrichment(shorter_pathways$ICD, signed_P_values_P3)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.1,1.1,0.1), limits = c(-0.1, 1.0))+theme_bw()
fgsea_shorter3

# Figure 2C Here
volcanoPlotWithList(Pat3_FL_end_Pat1_FL_end, l2fc_cutoff = 1, fdr = 0.05, title="Pat3 end vs Pat1 end", gene_list = pso_genes)

# Figure 2 D
# Extract genes which are differentially expressed between the end of Patient1FL and end of Patietn3 FL
significant_pso_genes <- Pat3_FL_end_Pat1_FL_end %>% filter(padj < 0.05 & log2FoldChange > 1 & Gene_name %in% pso_genes) %>% dplyr::select(Gene_name)
dds_Aldara <- DESeqDataSetFromMatrix(countData = round(Aldara_raw_count_filt), colData = sample_combi, design = ~1)
dds_Aldara <- DESeq(dds_Aldara)
vst_Aldara <- varianceStabilizingTransformation(dds_Aldara)
vst_Aldara <- assay(vst_Aldara)
vst_Aldara <- data.frame(vst_Aldara)
vst_Aldara$ensembl_id <- rownames(vst_Aldara)
vst_Aldara <- merge(vst_Aldara, gene_info[c("ensembl_id", "Gene_name")], by="ensembl_id")
dim(vst_Aldara) # 17153 x 9

vst_Aldara_pso <- data.frame(vst_Aldara) %>% filter(Gene_name %in% significant_pso_genes$Gene_name)
rownames(vst_Aldara_pso) <- vst_Aldara_pso$Gene_name
vst_Aldara_pso$Patient1 <- vst_Aldara_pso$MUC20494 - vst_Aldara_pso$MUC20493
vst_Aldara_pso$Patient3 <- vst_Aldara_pso$MUC20490 - vst_Aldara_pso$MUC20488
vst_Aldara_pso_Heatmap <- vst_Aldara_pso %>% dplyr::select(c(Patient1, Patient3))
vst_Aldara_pso_Heatmap <- data.frame(vst_Aldara_pso_Heatmap)
vst_Aldara_pso_Heatmap <- vst_Aldara_pso_Heatmap[order(vst_Aldara_pso_Heatmap$Patient3, decreasing = T),,drop=F]
Heatmap(data.matrix(vst_Aldara_pso_Heatmap), col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), show_row_dend = F, show_column_dend = F, cluster_rows = F, row_names_gp = gpar(fontsize = 6),
        height = unit(10, "cm"), width = unit(5, "mm"),column_names_gp = gpar(fontsize=6), heatmap_legend_param = list(title="Difference between end and start of treatment [vst counts]"))

# Figure 2 E
selected_genes <- c("KRT6C", "LCE3C", "HRNR", "IL19", "CXCL8", "TNIP3", "KLK13", "TMPRSS11D", "NOS2", "TCN1")
vst_Aldara_time <- data.frame(vst_Aldara) %>% filter(Gene_name %in% selected_genes)
rownames(vst_Aldara_time) <- vst_Aldara_time$Gene_name
vst_Aldara_time <- vst_Aldara_time %>% dplyr::select(c(MUC20488, MUC20490, MUC20492))
time_data <- t(vst_Aldara_time)
time_data <- data.frame(time_data)
time_data$sample <- rownames(time_data)
time_data <- merge(time_data, sample_combi[,c("sample","time")], by="sample")
time_data <- time_data[c(1,3,2),]
for (gene in selected_genes){
p <- ggplot(data = time_data, aes(x=factor(time, levels=c("start", "mid", "end")), y=!!sym(gene), group=1))+geom_line(linetype="dashed")+geom_point()+theme_bw()+ggtitle("VST counts from Day 0 to Day 28")+ylim(0,20)
print(p)
}

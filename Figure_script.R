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
library(DESeq2)
library(biomaRt)
library(dplyr)
library(enrichplot)
setwd("/Users/martin.meinel/Documents/Publication Code/Aldara_Recall")

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
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    )+
    ggtitle(title)
}

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

# Aldara samples alone filtered
Aldara_raw_count_filt <- Aldara_raw_count %>% filter(ensembl_id %in% filtered_named_TPM_counts_BRAIN_Aldara$GeneID)
dim(Aldara_raw_count_filt) # 17153 x 8

# Brain and Aldara samples filtered
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

#end1_FL: 
end1_FL_Aldara <- aldara_samples %>% 
  dplyr::filter(t_area_p == "end1_FL") %>% 
  dplyr::select(sampleID, dataset)

#end3_FL: End of Former lesional Patient
end3_FL_Aldara <- aldara_samples %>% 
  dplyr::filter(t_area_p == "end3_FL") %>% 
  dplyr::select(sampleID, dataset)

end3_NL_Aldara <- aldara_samples %>%
  dplyr::filter(t_area_p == "end3_NL") %>%
  dplyr::select(sampleID, dataset)

aldara_brain_comparison <- function(brain_samples, aldara_sample, Brain_aldara_counts, genes, genedata=geneDescription){
  aldara_brain_samples <- rbind(brain_samples, aldara_sample)
  aldara_brain_samples <- aldara_brain_samples[order(aldara_brain_samples$sampleID),]
  rownames(aldara_brain_samples) <- aldara_brain_samples$sampleID
  
  sorted_aldara_brain_counts <- BRAIN_Aldara_counts[,colnames(BRAIN_Aldara_counts) %in% aldara_brain_samples$sampleID]
  sorted_aldara_brain_counts <- sorted_aldara_brain_counts[,order(colnames(sorted_aldara_brain_counts))]
  assert(colnames(sorted_aldara_brain_counts) == aldara_brain_samples$sampleID)
  # Create the DGE
  y <- DGEList(counts = sorted_aldara_brain_counts,
               genes = genes,
               samples = aldara_brain_samples)
  dataset <- aldara_brain_samples$dataset
  design <- model.matrix(~0 + dataset)
  keep <- filterByExpr(y, design = design)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design=design)
  fit <- glmFit(y, design=design)
  lrt <- glmLRT(fit, contrast = c(1, -1))
  result_all <- as.data.frame(topTags(lrt, n = 18000, adjust.method = "BH", sort.by = "PValue", p.value = 1)) %>% 
    dplyr::mutate(DGE = case_when(FDR < 0.05 & logFC > 1 ~ "up",
                                  FDR < 0.05 & logFC < 1 ~ "down",
                                  T ~ "no"))%>%
    dplyr::inner_join(genedata, by = c("ensembl_id" = "ensembl_gene_id"))
  result_all <- result_all %>% dplyr::rename(padj=FDR, log2FoldChange=logFC, Gene_name=external_gene_name)
  return(result = result_all)
}

computeEnrichments <- function(comparison, pathways, plt_title){
  assert_colnames(comparison, colnames = c("padj", "log2FoldChange", "Gene_name"), only_colnames = F, quiet = F)
  # Preparation
  signed_p_values <- -log10(comparison$padj) * sign(comparison$log2FoldChange)
  names(signed_p_values) <- comparison$Gene_name
  signed_p_values <- signed_p_values[order(signed_p_values, decreasing = T)]
  signed_p_values <- signed_p_values[!is.na(signed_p_values)]
  # Computation
  fgsea_result <- fgsea(pathways, signed_p_values)
  # Get the lower limits
  p <- plotEnrichment(pathways$PSO, signed_p_values)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.5,1.1,0.1), limits = c(-0.5, 1.0), labels = scales::label_number(accuracy = 0.1))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p)
  p2 <- plotEnrichment(pathways$ICD, signed_p_values)+theme(text=element_text(size=8))+scale_y_continuous(breaks=seq(-0.2,1.1,0.1), limits = c(-0.2, 1.0), labels = scales::label_number(accuracy = 0.1))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(p2)
  return(fgsea_result)
}

#### Short Aldara - Brain Comparisons
Pat1FL_Brain_end <- aldara_brain_comparison(brain_samples = non_lesional_samples, aldara_sample = end1_FL_Aldara, 
                                            Brain_aldara_counts = Brain_aldara_counts, genes = genes_filtered_named_counts_Aldara_BRAIN)
Pat3NL_Brain_end <- aldara_brain_comparison(brain_samples = non_lesional_samples, aldara_sample = end3_NL_Aldara, 
                                            Brain_aldara_counts = Brain_aldara_counts, genes = genes_filtered_named_counts_Aldara_BRAIN) 
Pat3FL_Brain_end <- aldara_brain_comparison(brain_samples = non_lesional_samples, aldara_sample = end3_FL_Aldara, 
                                            Brain_aldara_counts = Brain_aldara_counts, genes = genes_filtered_named_counts_Aldara_BRAIN) 


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
my.contrasts <- makeContrasts(P3end_P1end = TPAend3_FL - TPAend1_FL, P3endFL_P3_end_NL = TPAend3_FL - TPAend3_NL,levels = design)
 # Different comparisons
Pat3_FL_end_Pat1_FL_end <- glmLRT(fit, contrast = my.contrasts[, "P3end_P1end"])
Pat3_FL_end_Pat1_FL_end <- topTags(Pat3_FL_end_Pat1_FL_end, n=length(Aldara_raw_count_filt$MUC20488),adjust.method = "BH")$table

Pat3_FL_end_Pat3_NL_end <- glmLRT(fit, contrast = my.contrasts[, "P3endFL_P3_end_NL"])
Pat3_FL_end_Pat3_NL_end <- topTags(Pat3_FL_end_Pat3_NL_end, n=length(Aldara_raw_count_filt$MUC20488), adjust.method = "BH")$table

Pat3_FL_end_Pat1_FL_end <- data.frame(merge(Pat3_FL_end_Pat1_FL_end, gene_info[,c("Gene_name", "ensembl_id")], by.x="genes", by.y="ensembl_id"))%>% dplyr::rename(log2FoldChange=logFC, padj=FDR)
Pat3_FL_end_Pat3_NL_end <- data.frame(merge(Pat3_FL_end_Pat3_NL_end, gene_info[,c("Gene_name", "ensembl_id")], by.x="genes", by.y="ensembl_id")) %>% dplyr::rename(log2FoldChange=logFC, padj=FDR)

# --------------------- ENRICHMENTS HERE --------------------------------------
# LOAD GENE SETS FOR PSORIASIS, ICD AND TRMS
pso_genes <- read_xlsx("./SciTranslMed_2014_Quaranta_Pso genes.xlsx")
pso_genes <- pso_genes %>% dplyr::select(gene)
pso_genes <- pso_genes[!is.na(pso_genes)]
pso_genes <- pso_genes[pso_genes %in% filtered_named_TPM_counts_BRAIN_Aldara$external_gene_name] # 81
ICD_genes <- read_xlsx("./Aldara_gene lists.xlsx", sheet="PL top genes") # 1072
# filter for protein coding genes
ICD_genes_merged <- ICD_genes %>% filter(hgnc_symbol %in% filtered_named_TPM_counts_BRAIN_Aldara$external_gene_name) # 956
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

### All enrichments
# Figure 1D and Supplementary Figure 1B
patient1_end_enrichment <- computeEnrichments(Pat1FL_Brain_end, shorter_pathways, "Pat1FL_end")
patient3_end_enrichmentNL <- computeEnrichments(Pat3NL_Brain_end, shorter_pathways, "Pat3NL_end")
patient3_end_enrichmentFL <- computeEnrichments(Pat3FL_Brain_end, shorter_pathways, "Pat3FL_end")

# Supllementary Figure 1c here
volcanoPlotWithList(Pat3_FL_end_Pat1_FL_end, l2fc_cutoff = 1, fdr = 0.05, title="Pat3 end vs Pat1 end", gene_list = pso_genes)
volcanoPlotWithList(Pat3_FL_end_Pat3_NL_end, l2fc_cutoff = 1, fdr = 0.05, title = "Pat3 end FL vs. Pat3 end NL", gene_list = pso_genes)
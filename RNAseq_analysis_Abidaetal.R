#Abida et al. 2019 article.
#some of the patients present in Abida et al.2019 are also part of Robinson et al.2015

################ load libraries ################ 
suppressMessages(library(edgeR))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(remotes))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(sva))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(biomaRt))
suppressMessages(library(tibble))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(KEGGREST))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(msigdbr))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))

setwd("/Users/maurizio.aurora/HSR Global Dropbox/WORKSPACE/Bellone/Revision_CR/Chinnaiyan/")


# read Abida et al. 2019 counts file provided by A. Chinnaiyan & M. Cieslik
my_filecount = 'rna_tcounts_SU2C.csv'
# read Abida et al. 2019 supplementary data S01 metadata file
metadata = read.xlsx("pnas.1902651116.sd01.xlsx", "Sheet1")
rownames(metadata) <- metadata$Sample.ID
table(metadata$Pathology.Classification)
# keep only the lines associated with the tumor types of interest
# format metadata to allow better string matching to link samples in metadata and counts file
metadata <- metadata[metadata$Pathology.Classification %in% c("Adenocarcinoma", "Adenocarcinoma with NE features", "Small cell"), ]
# clean the metadata
metadata$Sample.ID <- gsub("^PRAD-", "", metadata$Sample.ID)
table(metadata$Pathology.Classification)
# remove from metadata rows where Subject ID is highlighted in red
exclude_ids <- c("1115016", "1115082", "1115083", "1115087",
                 "6115392", "6115400", "9115050", "9115059",
                 "9115065", "9115174", "9115108", "9115341")
metadata <- metadata[!metadata$Subject.ID %in% exclude_ids, ]
# format string
metadata$Sample.IDclean <- sub("-.*$", "", metadata$Sample.ID)
# inspect file
table(metadata$Pathology.Classification)
#Adenocarcinoma Adenocarcinoma with NE features 
#308                              14 
#Small cell 
#27 

# read counts
fCounts <- read.delim(file=my_filecount, header=TRUE,
                      check.names = F, sep = ",", row.names = 1)

# Combine the two vectors of IDs to link metadata to counts iDs
# serching for either Sample.IDclean or Subject.ID
all_ids <- unique(c(metadata$Sample.IDclean, metadata$Subject.ID))

# Create a pattern with all IDs joined by |
pattern <- paste(all_ids, collapse = "|")

# Keep rows where rownames contain any of the IDs
fCounts_filtered <- fCounts[str_detect(rownames(fCounts), pattern), ]
fCounts_filtered <- unique(fCounts_filtered)
# Check
head(fCounts_filtered[, 1:3])
nrow(fCounts_filtered)

# Add a column with the matched ID
matched_ids <- str_extract(rownames(fCounts_filtered), pattern)
fCounts_filtered$MatchedID <- matched_ids

# Checks
nrow(fCounts_filtered)
nrow(unique(fCounts_filtered))
ncol(fCounts_filtered)
tail(fCounts_filtered$MatchedID)
tail(fCounts_filtered [, 33095:33098])

# Initialize the new column in counts
#fCounts_filtered$Sample.ID <- NA
# Initialize the new column in counts
#fCounts_filtered$Subject.ID <- NA

# Case 1: MatchedID is already a Sample.ID
#fCounts_filtered$Sample.ID[fCounts_filtered$MatchedID %in% metadata$Sample.ID] <- 
#  fCounts_filtered$MatchedID[fCounts_filtered$MatchedID %in% metadata$Sample.ID]

# Case 2: MatchedID is a Sample.IDclean, map to corresponding Sample.ID
#sample_mask <- fCounts_filtered$MatchedID %in% metadata$Sample.IDclean
#fCounts_filtered$Sample.ID[sample_mask] <- metadata$Sample.ID[match(
#  fCounts_filtered$MatchedID[sample_mask], 
#  metadata$Sample.IDclean
#)]


# Case 1: MatchedID is already a Subject.ID
#fCounts_filtered$Subject.ID[fCounts_filtered$MatchedID %in% metadata$Subject.ID] <- 
#  fCounts_filtered$MatchedID[fCounts_filtered$MatchedID %in% metadata$Subject.ID]

# Case 2: MatchedID is a Sample.IDclean, map to corresponding Subject.ID
#sample_mask <- fCounts_filtered$MatchedID %in% metadata$Sample.IDclean
#fCounts_filtered$Subject.ID[sample_mask] <- metadata$Subject.ID[match(
#  fCounts_filtered$MatchedID[sample_mask], 
#  metadata$Sample.IDclean
#)]



# Initialize new columns
fCounts_filtered$Sample.ID <- NA
fCounts_filtered$Subject.ID <- NA


# Mask: which MatchedIDs are contained in metadata$Sample.ID
mask <- sapply(fCounts_filtered$MatchedID, function(id)
  any(grepl(id, metadata$Sample.ID))
)

# Fill both Sample.ID and Subject.ID for matched rows
matched_data <- t(sapply(fCounts_filtered$MatchedID[mask], function(id) {
  row <- metadata[grepl(id, metadata$Sample.ID), c("Sample.ID", "Subject.ID")][1, ]
  c(Sample.ID = row$Sample.ID, Subject.ID = row$Subject.ID)
}))

fCounts_filtered$Sample.ID[mask] <- matched_data[, "Sample.ID"]
fCounts_filtered$Subject.ID[mask] <- matched_data[, "Subject.ID"]


# -----------------------------
# Fallback: if Subject.ID is still NA but Sample.ID exists, try to fill Subject.ID
# Fallback: match MatchedID to metadata$Subject.ID for rows where Subject.ID is still NA
fallback_mask <- is.na(fCounts_filtered$Subject.ID) & !is.na(fCounts_filtered$MatchedID)

# Only proceed if there are any such rows
if (any(fallback_mask)) {
  fallback_data <- t(sapply(fCounts_filtered$MatchedID[fallback_mask], function(id) {
    row <- metadata[grepl(id, metadata$Subject.ID), c("Sample.ID", "Subject.ID")][1, ]
    c(Sample.ID = row$Sample.ID, Subject.ID = row$Subject.ID)
  }))
  
  fCounts_filtered$Sample.ID[fallback_mask] <- fallback_data[, "Sample.ID"]
  fCounts_filtered$Subject.ID[fallback_mask] <- fallback_data[, "Subject.ID"]
}


fCounts_filtered$Batch <- ifelse(
  fCounts_filtered$MatchedID %in% metadata$Subject.ID,
  "Batch1",
  "Batch2"
)



# create new metadata file wherealso complicated counts samples ID appear
subset <- fCounts_filtered[, c("MatchedID", "Subject.ID", "Sample.ID", "Batch")]
subset$ExtName <- rownames(subset)
subset_full <- merge(subset, metadata, by = "Subject.ID", all.x = TRUE)

#add batch as identified from PCA and counts sampleIDs (here called "ExtName")
#subset_full$Batch <- ifelse(is.na(subset_full$Sample.ID.x), "Batch1", "Batch2")

table(subset_full$Pathology.Classification)

subset_full$Condition <- ifelse(
  subset_full$Pathology.Classification == "Adenocarcinoma", "Adenocarcinoma",
  ifelse(subset_full$Pathology.Classification %in% c("Adenocarcinoma with NE features", "Small cell"), "NE", 
         subset_full$Pathology.Classification)
)


# clean counts from extra columns
fCounts_filtered$Subject.ID <- NULL
fCounts_filtered$MatchedID <- NULL
fCounts_filtered$Sample.ID <- NULL
fCounts_filtered$Batch <- NULL
# transpose counts table
fCounts_filtered_t <- t(fCounts_filtered)

# Connect to Ensembl
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")  # human

gene_ids <- rownames(fCounts_filtered_t)

# Get annotation
annotation_cols <- c('ensembl_gene_id', 'chromosome_name', 'start_position', 
                     'end_position', 'strand', 'gene_biotype', 'external_gene_name')

gene_ids <- colnames(fCounts_filtered)

annotation_data <- getBM(
  attributes = annotation_cols,
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)


# Set rownames as Ensembl IDs
rownames(fCounts_filtered_t) <- colnames(fCounts_filtered)

# Now columns are samples
colnames(fCounts_filtered_t) <- rownames(fCounts_filtered)

annotation_data$length <- annotation_data$end_position - annotation_data$start_position + 1
rownames(annotation_data) <- annotation_data$ensembl_gene_id

length(unique(annotation_data$ensembl_gene_id))
length(unique(annotation_data$external_gene_name))

# keep only Ensembl IDs in annotation_data
fCounts_filtered_t <- fCounts_filtered_t[rownames(fCounts_filtered_t) %in% annotation_data$ensembl_gene_id, ]

# Reorder annotation to match counts
annotation_data <- annotation_data[rownames(fCounts_filtered_t), ]

# Bind annotation columns to counts
fCounts_annotated <- cbind(annotation_data, fCounts_filtered_t)
#head(fCounts_annotated[fCounts_annotated$external_gene_name == "A4GNT", ])

# Check result
head(fCounts_annotated[, 1:12])

# Summarise counts
counts_sum <- fCounts_annotated %>%
  filter(!is.na(external_gene_name) & external_gene_name != "") %>%
  group_by(external_gene_name) %>%
  summarise(across(where(is.numeric), sum))

nrow(counts_sum)
# Keep metadata (first occurrence)
metadata <- fCounts_annotated %>%
  filter(!is.na(external_gene_name) & external_gene_name != "") %>%
  group_by(external_gene_name) %>%
  slice(1) %>%
  dplyr::select(external_gene_name, ensembl_gene_id, chromosome_name, gene_biotype)


# Merge counts and metadata
fCounts_annotated_clean <- left_join(metadata, counts_sum, by = "external_gene_name") %>%
  column_to_rownames("external_gene_name")
nrow(fCounts_annotated_clean)

head(colnames(fCounts_annotated_clean), 10)
head(fCounts_annotated_clean[rownames(fCounts_annotated_clean) == "A4GNT", ])

fCounts_annotated_clean <- fCounts_annotated_clean %>%
  rownames_to_column(var = "external_gene_name") %>%
  relocate(external_gene_name, .before = 1)

rownames(fCounts_annotated_clean) <- fCounts_annotated_clean$external_gene_name
annotation <- c('external_gene_name','ensembl_gene_id','chromosome_name','gene_biotype','start_position','end_position','strand', 'length')


fCountsData <- fCounts_annotated_clean[
  , 
  -which(
    tolower(names(fCounts_annotated_clean))
    %in% 
      tolower(annotation))]


fCountsAnnotation <- fCounts_annotated_clean[
  , 
  which(
    tolower(names(fCounts_annotated_clean))
    %in% 
      tolower(annotation))]


geneidColname <- 'external_gene_name'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCounts_annotated_clean) <- fCounts_annotated_clean[[geneidIdx]]

table(subset_full$Condition)

table(subset_full$Batch, subset_full$Condition)

exp_batch=as.factor(subset_full$Batch)

# correct batch effect with ComBatSeq
# Using DEseq2 covariates in exp design formula is not enough.
# convert to numeric matrix
#fCounts_annotated_cs <- ComBat_seq(data.matrix(fCountsData), batch=exp_batch) #combined

#sum(duplicated(rownames(fCounts_annotated_cs)))
# Or see which genes are duplicated
#rownames(fCounts_annotated_cs)[duplicated(rownames(fCounts_annotated_cs))]

#duplicated_rows <- fCounts_annotated_cs[duplicated(fCounts_annotated_cs), ]

# Sum across each row (gene)
#row_sums <- rowSums(duplicated_rows)

dir.create("/Users/maurizio.aurora/HSR Global Dropbox/WORKSPACE/Bellone/Revision_CR/Chinnaiyan/ok/Adenocarcinoma_with_NE_and_small_cell_nocorr", recursive = TRUE)
getwd()
setwd("/Users/maurizio.aurora/HSR Global Dropbox/WORKSPACE/Bellone/Revision_CR/Chinnaiyan/ok/Adenocarcinoma_with_NE_and_small_cell_nocorr")

# Total sum across all duplicated rows
#total_sum <- sum(row_sums)
#length(unique(rownames(duplicated_rows)))
fCounts_annotated_cs <- fCountsData
# ~ Batch + Condition
#nrow(unique(fCounts_annotated_cs))
#nrow(fCounts_annotated_cs)
SAVE_variable <- list()

variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM', 'expGenes_RPKM', 'allGenes_LogRPKM')

rownames(subset_full) <- subset_full$ExtName
metadata_d = subset_full

#Adenocarcinoma Adenocarcinoma with NE features                      Adenocarcinoma with NE features 
#254  9                              13 

fCountsData_d = fCounts_annotated_cs[,row.names(metadata_d)]

# all_counts
#fCountsData_d <- as.data.frame(apply(fCountsData_d, 2, function(x) as.numeric(trimws(x))))
#rownames(fCountsData_d) <- rownames(fCountsData_d[,row.names(metadata_d)])

y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)

SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

Nreplica = 22
# expGENES_counts
keep <- rowSums(cpm(y)>1)>=Nreplica
table(keep)
yf <- y[keep,]

# CPM
SAVE_variable[[variable2save_names[2]]] <- as.data.frame(yf$counts)

# log CPM
SAVE_variable[[variable2save_names[3]]] <- as.data.frame(cpm(yf, log=T))

#RPKM log
SAVE_variable[[variable2save_names[4]]] <- as.data.frame(rpkm(yf, log=T, gene.length =yf$genes$length))

#RPKM not log
SAVE_variable[[variable2save_names[5]]] <- as.data.frame(rpkm(yf, log=F, gene.length =yf$genes$length))

#RPKM log all genes
SAVE_variable[[variable2save_names[6]]] <- as.data.frame(rpkm(y, log=T, gene.length =y$genes$length))

filename_xls <- "COUNTS_Abida_2019.xlsx"

write.xlsx(SAVE_variable,
           file = filename_xls, 
           rowNames = T,
           asTable = T, 
           sheetName =variable2save_names)

# calculate 500 most variant genes  

y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$length)
keep <- rowSums(cpm(y)>1)>=Nreplica
yf <- y[keep,]
nrow(yf)
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
provaA <- names(vary)

TOP_N <- names(vary_s[1:N])
nrow(TOP_N)
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]
#PCA parameters
pcx = 1
pcy = 2
centering = TRUE
scaling = TRUE
# PCA
pca = prcomp(t(fCountsRPKMTOP), center=centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)

# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)

# apparently there is a batch effect
#metadata_d$Batch <- ifelse(is.na(metadata_d$Sample.ID.x), "Batch1", "Batch2")
# let's verify
score$ExtName = metadata_d$ExtName
score$Pathology.Classification = metadata_d$Pathology.Classification
score$batch <- metadata_d$Batch
score$Sample.Type <- metadata_d$Sample.Type


pca <- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                         color=Pathology.Classification, shape=batch)) +
  geom_point(size=7) +
  labs(x=xlab, y=ylab, title=paste("PC", pcx, " vs PC", pcy, " scoreplot", sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey") +
  geom_vline(xintercept=0, linetype="dashed", color="darkgrey") +
  scale_color_manual(values=c("Adenocarcinoma" = "pink", "Adenocarcinoma with NE features" = "#66B2FF", "Small cell" = "darkgreen")) +  # Manually setting colors for Condition
  theme(plot.title=element_text(color="black", size=26, face="bold.italic"),
        axis.text.x=element_text(angle=0, face="bold", color="black", size=22, hjust=.5),
        axis.title.x=element_text(face="bold", color="black", size=24),
        axis.text.y=element_text(angle=0, face="bold", color="black", size=22),
        axis.title.y=element_text(face="bold", color="black", size=24),
        legend.text=element_text(face="bold", color="black", size=18),
        legend.position="right",
        panel.background=element_rect(fill="white", colour="black", size=1, linetype="solid"))

pdf('PCA_top500rpkm.pdf', width = 14, height = 6)
pca
dev.off()


annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 2
cols = gg_color_hue(n)

dev.new(width = 4, height = 4)

mycolors_b <- c('red','blue'); names(mycolors_b) = levels(annotation_column$Condition)

ann_colors = list()
ann_colors = list(
  Condition = mycolors_b
)

options(repr.plot.width=12, repr.plot.height=10)

annotation_column <- metadata_d[, c("Pathology.Classification", "Batch", "Abiraterone.(ABI).and.Enzalutamide.(ENZA).Exposure.Status","Taxane.exposure.status"), drop = FALSE]

ann_colors = list(
  Pathology.Classification = c("Adenocarcinoma" = "pink", "Adenocarcinoma with NE features" = "lightblue", "Small cell" = "darkgreen"),
  Taxane.exposure.status = c("Naïve" = "red", "Exposed"="purple", "UNK"="grey" ),
  `Abiraterone.(ABI).and.Enzalutamide.(ENZA).Exposure.Status` = c("Naive" = "red", "Exposed"="purple", "UNK"="grey", "On treatment" ="blue"),
  Batch = c("Batch1" = "orange", "Batch2"="#FFFF66"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)

ncol(fCountsRPKMTOP)
HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = T,
                         show_colnames = T,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_500rpkm___.pdf',
                         width = 40, height = 20)


HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = F,
                         show_colnames = F,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_500rpkm__small.pdf',
                         width = 10, height = 10)



# perform DGE

f = "Condition"
seqc_pvalue = 0.01
dgeResults = list()

metadata_s <- metadata_d
metadata_s$Condition <- factor(metadata_s$Condition)

comparison = list(Condition = c("NE", "Adenocarcinoma"))

for (i in names(comparison)) {
  print(i)
  fCountsData_s = fCountsData_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_s,
    colData  = metadata_s,
    design   = as.formula('~Condition'))
  
  filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
  table(filter)
  ddsFiltered <- dds[filter,]
  dga <- DESeq(
    object = ddsFiltered,
    test = "Wald",
    fitType = "parametric",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf)
  
  alpha = 0.05
  print(paste(comparison[[i]][1],"_vs_",comparison[[i]][2],sep=''))
  dgeResults.tmp <- results(dga,
                            contrast             = c(f,comparison[[i]][1],comparison[[i]][2]),
                            cooksCutoff          = Inf,
                            independentFiltering = TRUE,
                            alpha                = alpha,
                            pAdjustMethod        = "BH")
  summary(dgeResults.tmp)
  dgeResults[[i]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),]
  
  #PCA
  vsd <- vst(dga, blind=FALSE)
  main.factor = "Condition"
  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=sex))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=Condition, label = name),
                     #geom_label_repel(data= pcaData, aes(PC1, PC2, color=Condition, label = str_sub(name,-5,-1)),
                     size = 6,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                     segment.color = 'grey50') +
    geom_point(size=6) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(paste("PCA",i)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_color_manual(values = c('salmon','grey'))
  pdf(paste('pca_AM',i,'.pdf',sep=''),width=8, height=6)
  print(PCA)
  dev.off()
  
  seqcUP = row.names(dgeResults[[i]])[dgeResults[[i]]$pvalue <= seqc_pvalue &
                                        !is.na(dgeResults[[i]]$padj)&
                                        dgeResults[[i]]$log2FoldChange > 1]
  
  if (length(seqcUP) > 0) {
    # print heatmap
    annotation_column <- metadata_s
    row.names(annotation_column) <- row.names(metadata_s)
    options(repr.plot.width=12, repr.plot.height=10)
    crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
    colors = crp(255)
    head(cpm(counts(dga)))
    HP <- pheatmap::pheatmap(cpm(counts(dga))[seqcUP,],
                             scale = 'row',
                             cluster_rows = T,
                             cluster_cols = T,
                             show_rownames = F,
                             cutree_cols = 2,
                             fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                             display_numbers = F,
                             col=colors,
                             filename = paste('Heatmap_seqcUP_',i,'.pdf'),
                             width = 10, height = 11 )
  }
  
}

#LFC > 0 (up)       : 3825, 22%
#LFC < 0 (down)     : 3097, 18%

#write.xlsx(as.data.frame(dgeResults), 'DGE_results.xlsx', rowNames = T)

f = 'DGE_results'

dir.create(f, showWarnings=TRUE, recursive=TRUE)

lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))

dgeResults_table = list()
dgeResults_table = lapply(
  names(dgeResults),
  function(x)
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = paste(f,'/DGE_results.xlsx', sep=''),
           rowNames = F,
           asTable = T,
           startRow = 1,
           sheetName = str_sub(names(dgeResults),1,31))


DEGs <- as.data.frame(dgeResults_table$Condition)
rownames(DEGs) <- DEGs$Geneid
colnames(DEGs) <- c("GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj")
df <- DEGs
rownames(df) <- df$GeneID


library(ggrepel)
# Set thresholds
logFC_threshold <- 1
pval_threshold <- 0.05

# Create -log10(padj) for the plot
df$neg_log10_padj <- -log10(df$padj)

# Define color based on thresholds
df$color <- "grey"
df$color[df$log2FoldChange > logFC_threshold & df$padj < pval_threshold] <- "red"
df$color[df$log2FoldChange < -logFC_threshold & df$padj < pval_threshold] <- "blue"

table(df$color)

pdf("Adenocarcinoma_with_NE_vs_Adenocarcinoma_volcano_l2fc1_labels_Abida_2019.pdf", )
# Create the volcano plot
ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.75) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(abs(log2FoldChange) > logFC_threshold & padj < pval_threshold, GeneID, "")), 
                  box.padding = 0.3, max.overlaps = 10)+
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")
dev.off()

getwd()


significant_genes <- df[df$padj < 0.05, ]
significant <- significant_genes[significant_genes$log2FoldChange > 1, ]
nrow(significant) #2063

significant_genes <- df[df$padj < 0.05, ]
significant_genes<- significant_genes[significant_genes$log2FoldChange < -1, ]
nrow(significant_genes) #1323


significant_genes <- df[df$padj < 0.05, ]
significant_up <- significant_genes[significant_genes$log2FoldChange > 0, ]
nrow(significant_up) #3825

significant_genes <- df[df$padj < 0.05, ]
significant_dw <- significant_genes[significant_genes$log2FoldChange < 0, ]
nrow(significant_dw) #3097


significant <- rbind(significant_up, significant_dw)
head(significant)
nrow(significant) #6922
significant_genes <- significant$GeneID

ccc <- fCountsRPKM[rownames(fCountsRPKM) %in% significant_genes,]

annotation_column_ <- metadata_d[, c("Condition","Batch")]

annotation_column <- metadata_d[, c("Pathology.Classification", "Batch", "Abiraterone.(ABI).and.Enzalutamide.(ENZA).Exposure.Status","Taxane.exposure.status", "Sample.Type", "Condition"), drop = FALSE]

unique(metadata_d$Taxane.exposure.status)


ann_colors = list(
  Sample.Type=c("Liver" = "red", "Lung" = "#66B2FF", "Prostate" = "darkgreen", "Adrenal" = "pink", "Bone"= "yellow", "LN" = "purple", "Other Soft tissue" = "orange"),
  Pathology.Classification = c("Adenocarcinoma" = "pink",  "Adenocarcinoma with NE features" = "66B2FF", "Small cell" = "darkgreen"),
  Taxane.exposure.status = c("Naïve" = "red", "Exposed"="purple", "UNK"="grey" ),
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"),
  `Abiraterone.(ABI).and.Enzalutamide.(ENZA).Exposure.Status` = c("Naive" = "red", "Exposed"="purple", "UNK"="grey", "On treatment" ="blue"),
  Batch = c("Batch1" = "orange", "Batch2"="#FFFF66"))

ann_colors = list(
  Batch = c("Batch1" = "yellow",
            "Batch2" = "brown"),
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"))



HP <- pheatmap::pheatmap(ccc,
                         scale = 'row',
                         annotation_col = annotation_column_,
                         cluster_rows = T, 
                         cluster_cols = T, 
                         annotation_colors = ann_colors, 
                         show_rownames=T,
                         show_colnames=T,
                         main = "DEGs",
                         #treeheight_col  = 0,
                         treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdYlBu")))(100),
                         filename = 'DEGs_heatmap_abida_2019.pdf',
                         width = 7, height = 350)

HP <- pheatmap::pheatmap(ccc,
                         scale = 'row',
                         annotation_col = annotation_column,
                         cluster_rows = T, 
                         cluster_cols = T, 
                         annotation_colors = ann_colors, 
                         show_rownames=F,
                         show_colnames=F,
                         main = "DEGs",
                         #treeheight_col  = 0,
                         treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colors,
                         filename = 'DEGs_heatmap_Abida_2019small.pdf')


# Kegg enrichment with cluster Profiler

df <- DEGs
rownames(df) <- df$GeneID

significant_genes <- df[df$padj < 0.05, ]
significant_up <- significant_genes[significant_genes$log2FoldChange > 0, ]
up <- rownames(significant_up)
length(up) #3825


significant_genes <- df[df$padj < 0.05, ]
significant_dw <- significant_genes[significant_genes$log2FoldChange < 0, ]
dw <- rownames(significant_dw)
length(dw) #3097


upregulated_genes_entrez <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg_enrich_up <- enrichKEGG(gene         = upregulated_genes_entrez$ENTREZID,
                             organism     = 'hsa',        # human
                             pvalueCutoff = 0.05)

# View pathways
kegg_enrich_up$Description

dwregulated_genes_entrez <- bitr(dw, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg_enrich_dw <- enrichKEGG(gene         = dwregulated_genes_entrez$ENTREZID,
                             organism     = 'hsa',        # human
                             pvalueCutoff = 0.05)


go_bp_enrich_up <- enrichGO(
  gene          = upregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)

go_bp_enrich_dw <- enrichGO(
  gene          = dwregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)


go_cc_enrich_up <- enrichGO(
  gene          = upregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)

go_cc_enrich_dw <- enrichGO(
  gene          = dwregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)


# Prepare data
df_plot <- as.data.frame(kegg_enrich_dw) %>%
  mutate(logp = -log10(p.adjust)) %>%
  arrange(logp) %>%
  mutate(Description = fct_reorder(Description, logp))

write.xlsx(df_plot, "Dw_KEGG_in_Adenocarcinoma_with_NE_and_small_cell_vs_Adeno.xlsx")


pdf("Dw_KEGG_in_Adenocarcinoma_with_NE_and_small_cell_vs_Adeno.pdf", 10 ,10)
ggplot(df_plot, aes(x = logp, y = Description, fill = logp)) +
  geom_bar(stat = "identity", color = "black", size = 0.6, width = 0.6) +  # black border
  scale_fill_gradient(
    low = "white",  # pale red
    high = "#DC143C"  # deep red
  ) +
  scale_x_continuous(expand = c(0,0)) +   # bars start at y-axis (no gap)
  labs(
    title = "Downregulated KEGG pathways \
    in Adeno with NE features + small cell  vs Adeno",
    x = expression(-log[10](adjusted~p~value)),
    y = NULL,
    fill = expression(-log[10](p.adjust))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.ticks   = element_line(color = "black", size = 0.6),   # black ticks
    axis.line    = element_line(color = "black", size = 0.8),   # black axes
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
  )
dev.off()


kegg_enrich_up$Description
pdf("KEGG_Enrichment_Dotplot_DW_batch.pdf")
dotplot(kegg_enrich_dw, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot DW") +
  theme_bw(base_size = 14)
dev.off()


pdf("KEGG_Enrichment_Dotplot_UP_batch.pdf")
dotplot(kegg_enrich_up, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot UP") +
  theme_bw(base_size = 14)
dev.off()


pdf("GOBP_Enrichment_Dotplot_DW_batch.pdf",10,10)
dotplot(go_bp_enrich_dw, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot UP") +
  theme_bw(base_size = 14)
dev.off()


pdf("GOBP_Enrichment_Dotplot_UP_batch.pdf",10,10)
dotplot(go_bp_enrich_up, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot UP") +
  theme_bw(base_size = 14)
dev.off()



# Get the genes in KEGG pathway hsa04512
kegg_e <- keggGet("hsa04512")[[1]]$GENE

# The returned vector has alternating EntrezID and gene symbol + description
# Extract the gene symbols
genes_entrez <- kegg_e[seq(1, length(kegg_e), 2)]
genes_info   <- kegg_e[seq(2, length(kegg_e), 2)]

# Parse out the gene symbols from the description strings
genes_symbols <- sapply(strsplit(genes_info, " "), `[`, 1)

# Create a data.frame
pathway_genes <- data.frame(entrez   = genes_entrez,
                            symbol   = genes_symbols,
                            stringsAsFactors = FALSE)

# Remove trailing semicolons
pathway_genes$symbol <- gsub(";$", "", pathway_genes$symbol)

# Check again
head(pathway_genes$symbol)
length(pathway_genes$symbol) #89

# keep u and down regulated DEGs
all_genes <- c(up, dw)
# for these get fCountsRPKMs
fCountsRPKM_DEGS <- fCountsRPKM[rownames(fCountsRPKM) %in% all_genes,]
head(fCountsRPKM_DEGS[, 1:3])

fCountsRPKM_DEGS_KEGG <- fCountsRPKM_DEGS[rownames(fCountsRPKM_DEGS) %in% pathway_genes$symbol,]

annotation_column__ <- as.data.frame(metadata_d[, c("Condition")])
rownames(annotation_column__) <- rownames(metadata_d)

ann_colors_ = list(
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
fCountsRPKM_DEGS_KEGG <- data.matrix(fCountsRPKM_DEGS_KEGG) 
rownames(fCountsRPKM_DEGS_KEGG)
colnames(annotation_column__) <- "Condition"

# Define your desired order
desired_order <- c("Adenocarcinoma", "NE")

# Reorder columns of the matrix based on that order
ordered_samples <- rownames(annotation_column__)[order(factor(annotation_column__$Condition,
                                                              levels = desired_order))]

# Reorder both matrix and annotation accordingly
fCountsRPKM_DEGS_KEGG <- fCountsRPKM_DEGS_KEGG[, ordered_samples]
annotation_column__ <- annotation_column__[ordered_samples, , drop = FALSE]


HP <- pheatmap::pheatmap(fCountsRPKM_DEGS_KEGG,
                         scale = 'row',
                         annotation_col = annotation_column__,
                         cluster_rows = T, 
                         cluster_cols = F, 
                         annotation_colors = ann_colors_, 
                         show_rownames=T,
                         show_colnames=F,
                         display_numbers = F, 
                         border_color=NA,
                         col=colors,
                         filename = 'DEGs_heatmap_ECM_KEGG_test_LFC_0_allDEGs_final.pdf',
                         width = 10, height = 8)


genes_GO0030198 <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0030198",
  keytype = "GOALL",    # all GO annotations
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
)


# Extract the gene symbols
genes_entrez <- genes_GO0030198$ENTREZID
genes_info   <- genes_GO0030198$SYMBOL

# Parse out the gene symbols from the description strings
genes_symbols <- unique(genes_info)
length(genes_symbols) #331

# keep u and down regulated DEGs
all_genes <- c(up, dw)
# for these get fCountsRPKMs
fCountsRPKM_DEGS <- fCountsRPKM[rownames(fCountsRPKM) %in% all_genes,]
head(fCountsRPKM_DEGS[, 1:3])

fCountsRPKM_GO0030198 <- fCountsRPKM_DEGS[rownames(fCountsRPKM_DEGS) %in% genes_info,]

annotation_column__ <- as.data.frame(metadata_d[, c("Condition")])
rownames(annotation_column__) <- rownames(metadata_d)

ann_colors_ = list(
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
fCountsRPKM_DEGS_KEGG <- data.matrix(fCountsRPKM_DEGS_KEGG) 
rownames(fCountsRPKM_DEGS_KEGG)
colnames(annotation_column__) <- "Condition"

# Define your desired order
desired_order <- c("Adenocarcinoma", "NE")

# Reorder columns of the matrix based on that order
ordered_samples <- rownames(annotation_column__)[order(factor(annotation_column__$Condition,
                                                              levels = desired_order))]

# Reorder both matrix and annotation accordingly
fCountsRPKM_GO0030198 <- fCountsRPKM_GO0030198[, ordered_samples]
annotation_column__ <- annotation_column__[ordered_samples, , drop = FALSE]


HP <- pheatmap::pheatmap(fCountsRPKM_GO0030198,
                         scale = 'row',
                         annotation_col = annotation_column__,
                         cluster_rows = T, 
                         cluster_cols = F, 
                         annotation_colors = ann_colors_, 
                         show_rownames=T,
                         show_colnames=F,
                         display_numbers = F, 
                         border_color=NA,
                         col=colors,
                         filename = 'DEGs_heatmap_ECM_GO0030198_test_LFC_0_allDEGs_final.pdf',
                         width = 10, height = 12)




split_gene_ids <- strsplit(go_bp_enrich_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Hs.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_dw <- as.data.frame(go_bp_enrich_dw)


go_enrichment_BP_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_dw, file = "go_enrichment_BP_dw_results.xlsx")


split_gene_ids <- strsplit(go_cc_enrich_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Hs.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  return(mapped_symbols)
}


# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_cc_enrich_dw <- as.data.frame(go_cc_enrich_dw)


go_cc_enrich_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_cc_enrich_dw, file = "go_enrichment_CC_dw_results.xlsx")



###############################################################################
###############################################################################
###############################################################################
# Analysis with LFC > 1 or < -1

dir.create("/Users/maurizio.aurora/HSR Global Dropbox/WORKSPACE/Bellone/Revision_CR/Chinnaiyan/ok/Adenocarcinoma_with_NE_and_small_cell_nocorr/LFC1", recursive = TRUE)
getwd()
setwd("/Users/maurizio.aurora/HSR Global Dropbox/WORKSPACE/Bellone/Revision_CR/Chinnaiyan/ok/Adenocarcinoma_with_NE_and_small_cell_nocorr/LFC1")

df <- DEGs
rownames(df) <- df$GeneID

# Kegg enrichment with cluster Profiler

significant_genes <- df[df$padj < 0.05, ]
significant_up <- significant_genes[significant_genes$log2FoldChange > 1, ]
up <- rownames(significant_up)
length(up) #2063

significant_genes <- df[df$padj < 0.05, ]
significant_dw <- significant_genes[significant_genes$log2FoldChange < -1, ]
dw <- rownames(significant_dw)
length(dw) #1323

upregulated_genes_entrez <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg_enrich_up <- enrichKEGG(gene         = upregulated_genes_entrez$ENTREZID,
                             organism     = 'hsa',        # human
                             pvalueCutoff = 0.05)

dwregulated_genes_entrez <- bitr(dw, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kegg_enrich_dw <- enrichKEGG(gene         = dwregulated_genes_entrez$ENTREZID,
                             organism     = 'hsa',        # human
                             pvalueCutoff = 0.05)


go_bp_enrich_up <- enrichGO(
  gene          = upregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)

go_bp_enrich_dw <- enrichGO(
  gene          = dwregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)


go_cc_enrich_up <- enrichGO(
  gene          = upregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)

go_cc_enrich_dw <- enrichGO(
  gene          = dwregulated_genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   # optional: converts ENTREZ → gene symbols
)



# View pathways
kegg_enrich_dw$Description

# Prepare data
df_plot <- as.data.frame(kegg_enrich_dw) %>%
  mutate(logp = -log10(p.adjust)) %>%
  arrange(logp) %>%
  mutate(Description = fct_reorder(Description, logp))

write.xlsx(df_plot, "Dw_KEGG_in_Adenocarcinoma_with_NE_and_small_cell_vs_Adeno.xlsx")
getwd()

pdf("Dw_KEGG_in_Adenocarcinoma_with_NE_and_small_cell_vs_Adeno.pdf", 10 ,10)
ggplot(df_plot, aes(x = logp, y = Description, fill = logp)) +
  geom_bar(stat = "identity", color = "black", size = 0.6, width = 0.6) +  # black border
  scale_fill_gradient(
    low = "white",  # pale red
    high = "#DC143C"  # deep red
  ) +
  scale_x_continuous(expand = c(0,0)) +   # bars start at y-axis (no gap)
  labs(
    title = "Downregulated KEGG pathways \
    in Adeno with NE features + Small cell vs Adeno",
    x = expression(-log[10](adjusted~p~value)),
    y = NULL,
    fill = expression(-log[10](p.adjust))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.ticks   = element_line(color = "black", size = 0.6),   # black ticks
    axis.line    = element_line(color = "black", size = 0.8),   # black axes
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
  )
dev.off()


kegg_enrich_up$Description
pdf("KEGG_Enrichment_Dotplot_DW.pdf")
dotplot(kegg_enrich_dw, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot DW") +
  theme_bw(base_size = 14)
dev.off()


pdf("KEGG_Enrichment_Dotplot_UP.pdf")
dotplot(kegg_enrich_up, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot UP") +
  theme_bw(base_size = 14)
dev.off()


pdf("GOBP_Enrichment_Dotplot_DW.pdf",10,10)
dotplot(go_bp_enrich_dw, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot DW") +
  theme_bw(base_size = 14)
dev.off()


pdf("GOBP_Enrichment_Dotplot_UP.pdf",10,10)
dotplot(go_bp_enrich_up, showCategory = 30) + 
  ggtitle("KEGG Enrichment Dotplot UP") +
  theme_bw(base_size = 14)
dev.off()

# Get the genes in KEGG pathway hsa04512
kegg_e <- keggGet("hsa04512")[[1]]$GENE

# The returned vector has alternating EntrezID and gene symbol + description
# Extract the gene symbols
genes_entrez <- kegg_e[seq(1, length(kegg_e), 2)]
genes_info   <- kegg_e[seq(2, length(kegg_e), 2)]

# Parse out the gene symbols from the description strings
genes_symbols <- sapply(strsplit(genes_info, " "), `[`, 1)

# Create a data.frame
pathway_genes <- data.frame(entrez   = genes_entrez,
                            symbol   = genes_symbols,
                            stringsAsFactors = FALSE)


# Remove trailing semicolons
pathway_genes$symbol <- gsub(";$", "", pathway_genes$symbol)

# Check again
head(pathway_genes$symbol)
length(pathway_genes$symbol) #89

# keep u and down regulated DEGs
all_genes <- c(up, dw)
# for these get fCountsRPKMs
fCountsRPKM_DEGS <- fCountsRPKM[rownames(fCountsRPKM) %in% all_genes,]
head(fCountsRPKM_DEGS[, 1:3])

fCountsRPKM_DEGS_KEGG <- fCountsRPKM_DEGS[rownames(fCountsRPKM_DEGS) %in% pathway_genes$symbol,]

annotation_column__ <- as.data.frame(metadata_d[, c("Condition")])
rownames(annotation_column__) <- rownames(metadata_d)

ann_colors_ = list(
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
fCountsRPKM_DEGS_KEGG <- data.matrix(fCountsRPKM_DEGS_KEGG) 
rownames(fCountsRPKM_DEGS_KEGG)
colnames(annotation_column__) <- "Condition"

# Define your desired order
desired_order <- c("Adenocarcinoma", "NE")

# Reorder columns of the matrix based on that order
ordered_samples <- rownames(annotation_column__)[order(factor(annotation_column__$Condition,
                                                              levels = desired_order))]

# Reorder both matrix and annotation accordingly
fCountsRPKM_DEGS_KEGG <- fCountsRPKM_DEGS_KEGG[, ordered_samples]
annotation_column__ <- annotation_column__[ordered_samples, , drop = FALSE]


HP <- pheatmap::pheatmap(fCountsRPKM_DEGS_KEGG,
                         scale = 'row',
                         annotation_col = annotation_column__,
                         cluster_rows = T, 
                         cluster_cols = F, 
                         annotation_colors = ann_colors_, 
                         show_rownames=T,
                         show_colnames=F,
                         display_numbers = F, 
                         border_color=NA,
                         col=colors,
                         filename = 'DEGs_heatmap_ECM_KEGG_test_LFC_1_allDEGs_final.pdf',
                         width = 10, height = 8)


genes_GO0030198 <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0030198",
  keytype = "GOALL",    # all GO annotations
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
)


# Extract the gene symbols
genes_entrez <- genes_GO0030198$ENTREZID
genes_info   <- genes_GO0030198$SYMBOL

# Parse out the gene symbols from the description strings
genes_symbols <- unique(genes_info)
length(genes_symbols) #331

# keep u and down regulated DEGs
all_genes <- c(up, dw)
# for these get fCountsRPKMs
fCountsRPKM_DEGS <- fCountsRPKM[rownames(fCountsRPKM) %in% all_genes,]
head(fCountsRPKM_DEGS[, 1:3])

fCountsRPKM_GO0030198 <- fCountsRPKM_DEGS[rownames(fCountsRPKM_DEGS) %in% genes_info,]

annotation_column__ <- as.data.frame(metadata_d[, c("Condition")])
rownames(annotation_column__) <- rownames(metadata_d)

ann_colors_ = list(
  Condition = c("Adenocarcinoma" = "purple", "NE"="darkgreen"))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
fCountsRPKM_DEGS_KEGG <- data.matrix(fCountsRPKM_DEGS_KEGG) 
rownames(fCountsRPKM_DEGS_KEGG)
colnames(annotation_column__) <- "Condition"

# Define your desired order
desired_order <- c("Adenocarcinoma", "NE")

# Reorder columns of the matrix based on that order
ordered_samples <- rownames(annotation_column__)[order(factor(annotation_column__$Condition,
                                                              levels = desired_order))]

# Reorder both matrix and annotation accordingly
fCountsRPKM_GO0030198 <- fCountsRPKM_GO0030198[, ordered_samples]
annotation_column__ <- annotation_column__[ordered_samples, , drop = FALSE]


HP <- pheatmap::pheatmap(fCountsRPKM_GO0030198,
                         scale = 'row',
                         annotation_col = annotation_column__,
                         cluster_rows = T, 
                         cluster_cols = F, 
                         annotation_colors = ann_colors_, 
                         show_rownames=T,
                         show_colnames=F,
                         #main = "DEGs",
                         #treeheight_col  = 0,
                         #treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colors,
                         filename = 'DEGs_heatmap_ECM_GO0030198_test_LFC_1_allDEGs_final.pdf',
                         width = 10, height = 12)

genes_order <- c("CREB3L1", "SLC2A10", "ADAMTS1", "ADAM15", "KLK4", "TPSAB1", "LCP1","FMOD",
    "ITGA8", "MYH11", "TNXB", "ZNF469", "ADAMTS9", "COMP", "AEBP1",
    "COL12A1","IL6", "MMP1", "ADAMTS8", "LAMA1", "MATN4", "KLK7", 
    "WNT3A", "MMP8", "CTSG", "ELANE", "PTX3", "VIT", "ITGB3", "DMP1",
    "IBSP", "ADAMTS3", "CYP1B1", "IHH","ADTRP", "AXIN2", "APLP1", "CARMIL2", "SMPD3", "WASHC1",
    "SMOC1", "TMPRSS6", "PLG", "VTN", "AGT",
    "SERPINF2",  "MMP16", "ADAMTS14", "ADAMTS18", "LARGE1", "HAS3",
    "COL4A6", "COL4A4", "COL6A6", "TNR")

idx <- match(genes_order, rownames(fCountsRPKM_GO0030198))
fCountsRPKM_ordered <- fCountsRPKM_GO0030198[idx[!is.na(idx)], , drop = FALSE]

HP <- pheatmap::pheatmap(
  fCountsRPKM_ordered,
  scale = "row",
  annotation_col = annotation_column__,
  cluster_rows = FALSE,   # IMPORTANT: keep as given order
  cluster_cols = FALSE,
  annotation_colors = ann_colors_,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0,
  display_numbers = FALSE,
  border_color = NA,
  col = colors,
  filename = "DEGs_heatmap_ECM_GO0030198_test_LFC_1_allDEGs_ordered.pdf",
  width = 10,
  height = 12)

idx <- match(genes_order, rownames(fCountsRPKM_GO0030198))
fCountsRPKM_ordered <- fCountsRPKM_GO0030198[idx[!is.na(idx)], , drop = FALSE]

HP <- pheatmap::pheatmap(
  fCountsRPKM_ordered,
  scale = "row",
  annotation_col = annotation_column__,
  cluster_rows = FALSE,   # IMPORTANT: keep as given order
  cluster_cols = FALSE,
  annotation_colors = ann_colors_,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0,
  display_numbers = FALSE,
  border_color = NA,
  col = colors,
  filename = "DEGs_heatmap_ECM_GO0030198_test_LFC_1_allDEGs_ordered.pdf",
  width = 10,
  height = 12)


idx <- match(genes_order, rownames(fCountsRPKM_GO0030198))
fCountsRPKM_ordered <- fCountsRPKM_GO0030198[idx[!is.na(idx)], , drop = FALSE]

HP <- pheatmap::pheatmap(
  fCountsRPKM_ordered,
  scale = "row",
  annotation_col = annotation_column__,
  cluster_rows = FALSE,   # IMPORTANT: keep as given order
  cluster_cols = FALSE,
  annotation_colors = ann_colors_,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0,
  display_numbers = FALSE,
  border_color = NA,
  col = colors,
  filename = "DEGs_heatmap_ECM_GO0030198_test_LFC_1_allDEGs_ordered.pdf",
  width = 10,
  height = 12
)


split_gene_ids <- strsplit(go_bp_enrich_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Hs.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  return(mapped_symbols)
}

# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_enrichment_BP_dw <- as.data.frame(go_bp_enrich_dw)


go_enrichment_BP_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_enrichment_BP_dw, file = "go_enrichment_BP_dw_results.xlsx")


split_gene_ids <- strsplit(go_cc_enrich_dw$geneID, "/")
map_gene_ids <- function(ids) {
  mapped_symbols <- mapIds(org.Hs.eg.db, 
                           keys = ids, 
                           column = "SYMBOL", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  return(mapped_symbols)
}


# Apply this function to each list of Entrez IDs
mapped_gene_symbols <- lapply(split_gene_ids, map_gene_ids)
go_cc_enrich_dw <- as.data.frame(go_cc_enrich_dw)


go_cc_enrich_dw$gene_symbol <- sapply(mapped_gene_symbols, function(symbols) {
  paste(symbols, collapse = "/")  # You can change '/' to ',' if preferred
})

write.xlsx(go_cc_enrich_dw, file = "go_enrichment_CC_dw_results.xlsx")

getwd()





# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$GeneID
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

names(gene_list)
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")

head(msigdbr_df)
df <- msigdbr_df[,c("gs_name","gene_symbol")]
msigdbr_df <- df 

head(msigdbr_df)

# Perform GSEA analysis 
gsea_results_H <- GSEA(gene_list, 
                       TERM2GENE = msigdbr_df,  # Replace 'C1' with the actual gene set data frame
                       verbose = FALSE, 
                       pvalueCutoff = 0.05, 
                       maxGSSize = 5000, 
                       minGSSize = 5, 
                       pAdjustMethod = "BH", 
                       eps = 0)



write.xlsx(gsea_results_H[gsea_results_H$Description], "H_gsea.xlsx")




msigdbr_df <- msigdbr(species = "Homo sapiens", category = "C2")

head(msigdbr_df)
df <- msigdbr_df[,c("gs_name","gene_symbol")]
msigdbr_df <- df 

head(msigdbr_df)

# Perform GSEA analysis (ensure TERM2GENE, e.g., 'C1', is properly defined)
gsea_results_C2 <- GSEA(gene_list, 
                        TERM2GENE = msigdbr_df,  # Replace 'C1' with the actual gene set data frame
                        verbose = FALSE, 
                        pvalueCutoff = 0.05, 
                        maxGSSize = 5000, 
                        minGSSize = 5, 
                        pAdjustMethod = "BH", 
                        eps = 0)



pdf('C2_KEGG_ECM_RECEPTOR_INTERACTION.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "KEGG_ECM_RECEPTOR_INTERACTION"), 
          title = "KEGG_ECM_RECEPTOR_INTERACTION")
dev.off()


pdf('C2_KEGG_FOCAL_ADHESION.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "KEGG_FOCAL_ADHESION"), 
          title = "KEGG_FOCAL_ADHESION")
dev.off()


pdf('C2_WP_FOCAL_ADHESION_PI3KAKTMTORSIGNALING.pdf', width = 6, height = 5)
gseaplot2(gsea_results_C2, geneSetID = which(gsea_results_C2$Description == "WP_FOCAL_ADHESION_PI3KAKTMTORSIGNALING"), 
          title = "WP_FOCAL_ADHESION_PI3KAKTMTORSIGNALING")
dev.off()


write.xlsx(gsea_results_C2[gsea_results_C2$Description], "C2_gsea.xlsx")

###############################################################################
###############################################################################




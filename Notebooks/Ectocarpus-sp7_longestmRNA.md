---
Master table of Ectocarpus sp7
---

This notebook aims to keep track of how I made the MASTER TABLE that is used for all downstream analyses.
*gff* used: `/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/annotation/reference/Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gff` + `/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/annotation/alternative_versions/Ec32_V5.final.bac_filtered.liniar_T2T_Ec32_v3.gt.merged_overlap.longest.mRNA_200bp.gff` = `Ec32_V5.longest_mRNA.SDR_F.mod.gff` custom version and *fasta*: `/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/genome/reference/Ec32_V5.final.bac_filtered.liniar_T2T.reference.fasta`.

#### Set up

```{r, message = FALSE}
# Loading lib -------------------------------------------------------------
library(rtracklayer)
library(dplyr)
library(dlookr)
library(tidyr)
library(stringr)
library(dlookr)
library(data.table)
library(GenomicRanges)
library(readr)
library(withr)

```


# ---------------------------------------------
# BLOCK1: Fill in gene_id gene_start gene_end with gff file 
# ---------------------------------------------

```{r}
## Extract gene_id, gene_start, and gene_end from gff file -------------------

# Define file path
gff_file <- "./data/gff_files/Ec32_V5.longest_mRNA.SDR_F.AGAT.mod.gtf"
# gff_file <- "./data/gff_files/Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gff"
# Import GFF data
gff_data <- readGFF(gff_file)

# Filter for mRNA features
gff_data$type <- as.character(gff_data$type)
mRNA <- gff_data[gff_data$type == "mRNA", ]

# Convert to GRanges object
mRNA_gr <- GRanges(
    seqnames = mRNA$seqid,
    ranges = IRanges(start = mRNA$start, end = mRNA$end),
    strand = mRNA$strand, 
    gene_id = mRNA$ID)


## Fill in the values in the master table -----------------------------------
master_table <- as.data.table(mRNA_gr) %>%  distinct()
colnames(master_table) <- c("seq_id", "gene_start", "gene_end", 
                            "width", "strand", "gene_id")


print(master_table)

```


# ---------------------------------------------
# BLOCK2: Fetch TPM from salmon data
# ---------------------------------------------

```{r}
## Step 1: Read the Salmon output --------------------------------------------
salmon_file <- "./data/my_RNAseq_data/Ectocarpus-sp7/salmon.merged.gene_tpm.tsv" 
salmon_data <- read.table(salmon_file, header = TRUE, sep = "\t")
colnames(salmon_data)[1] <- 'trad_id'
salmon_data <- salmon_data %>%
  select(-trad_id) %>%  # Remove unnecessary columns
  rename(gene_id = gene_name)  # Rename column
  # group_by(gene_id) %>%  # Group by gene_id
  # summarise(across(everything(), sum, na.rm = TRUE)) # Sum all numeric columns

## Step 2: Calculate the mean to create the TPM column ----
salmon_data$TPM.Female <- rowMeans(salmon_data[, 2:3], na.rm = TRUE)
salmon_data$TPM.Male <- rowMeans(salmon_data[, 4:5], na.rm = TRUE)

## Step 3: Calculate log2(TPM + 1) and fill in the respective column ---------------
salmon_data$log2_TPM_plus_1.Female <- log2(salmon_data$TPM.Female + 1)
salmon_data$log2_TPM_plus_1.Male <- log2(salmon_data$TPM.Male + 1)

salmon_data <- salmon_data %>% select(gene_id, TPM.Female, TPM.Male, log2_TPM_plus_1.Female, log2_TPM_plus_1.Male)
head(salmon_data)

# ## Step 4: Translate trad_id to gene_id thanks to a dictionnary (made with Python script not here) ----------------------- 
# dict_trad <- read_csv('./data/Ectocarpus-sp7-mRNA_dict2.csv', show_col_types = FALSE)
# colnames(dict_trad) <- c('chromid', 'start', 'end', 'gene_id','trad_id')
# head(dict_trad)

## Step 5: Merge TPM values with the master table based on trad_id -----------
master_table <- left_join(master_table, salmon_data, by = "gene_id")
head(master_table %>% select(seq_id, gene_id, TPM.Female, TPM.Male, log2_TPM_plus_1.Female, log2_TPM_plus_1.Male))
## Verify the merge was successful and there are no NA values where there shouldn't be ------
if (any(is.na(master_table$TPM.Female | master_table$TPM.Male))) {
  warning("Some geneIDs in the master table do not have matching TPM values in the Salmon output.")
  missing_gene_ids <- master_table$gene_id[!(master_table$gene_id %in% salmon_data$gene_id)]
print(missing_gene_ids)
print('I am removing them from the rest of the study. Only the first transcript is kept, as I made a custom gtf file I should have remove all the transcripts from sdr female that are not .1 .')
master_table <- master_table %>% filter(gene_id %in% salmon_data$gene_id)

}
```

# ---------------------------------------------
# BLOCK3: K79 peak annotation
# ---------------------------------------------
### FEMALE

```{r, warning=FALSE}
# Load ChIP-seq peak file
peaks <- read.table("./data/peak_K79/Ectocarpus-sp7_Ec25_FEMALE_H3K79me2_merged-peaks_annotatePeaks_gene.txt", header = TRUE, sep = '\t')
peaks_gr <- GRanges(seqnames = peaks$Chr,
                    ranges = IRanges(start = peaks$Start, 
                                     end = peaks$End
                                     ))

# Find overlaps between genes and peaks
overlaps <- findOverlaps(mRNA_gr, peaks_gr)

# Create a new column in the master table to indicate if a peak overlaps
master_table$peak_K79.Female <- ifelse(seq_along(master_table$gene_id) %in% queryHits(overlaps), TRUE, FALSE)

```

### MALE
```{r, warning=FALSE}
# Load ChIP-seq peak file
peaks <- read.table("./data/peak_K79/Ectocarpus-sp7_Ec32_MALE_H3K79me2_merged-peaks_annotatePeaks_gene.txt", header = TRUE, sep = '\t')
peaks_gr <- GRanges(seqnames = peaks$Chr,
                    ranges = IRanges(start = peaks$Start, 
                                     end = peaks$End
                                     ))

# Find overlaps between genes and peaks
overlaps <- findOverlaps(mRNA_gr, peaks_gr)

# Create a new column in the master table to indicate if a peak overlaps
master_table$peak_K79.Male <- ifelse(seq_along(master_table$gene_id) %in% queryHits(overlaps), TRUE, FALSE)

```

# ---------------------------------------------
# BLOCK4: Sex-bias
# ---------------------------------------------

```{r, warning=FALSE, message=FALSE}
## Step 1: Generates volcano.data - RNA-seq analysis of Ecto between male and female

#source("https://bioconductor.org/biocLite.R")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("dplyr")
library("tximport")
library("ComplexHeatmap")
library("ggplot2")
library("tidyr")

#### run DESeq2
## --- read in counts and info
countsToUse <- read.delim("./data/my_RNAseq_data/Ectocarpus-sp7/salmon.merged.gene_counts.tsv", header = T, as.is = T)
# countsToUse <- countsToUse[,c(1,2,3,4,6,5,7,8)] #if reorder needed
colnames(countsToUse)[c(3:6)] <- c("EcF1", "EcF2", "EcM1", "EcM2")
rownames(countsToUse) <-countsToUse$gene_id
countsToUse <- countsToUse[,-c(1,2)]
countsToUse <- as.matrix(countsToUse)
countsToUse <- subset(countsToUse, rowMax(countsToUse) >= 10)   #### this is optional

colData <- data.frame(genotype = c(rep("Female",2), rep("Male",2)), 
                      row.names = colnames(countsToUse) )

dds <- DESeqDataSetFromMatrix(round(countsToUse), colData = colData, design=~genotype)
dds <- estimateSizeFactors(dds)
run.dds <- DESeq(dds)

#pdf("/ebio/abt5_projects/", width = 8, height = 8)
rld <- rlog(run.dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$genotype <- factor(pcaData$genotype, c("Male","Female"))
plotPCA <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=3, aes(color=genotype)) +
  scale_color_manual(values =  c("#1f78b4", "#fb8072")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave(filename = "plotPCA_RNAseqx3_Ectocarpus-sp7-WT-lastassembly.pdf", 
       plot = plotPCA,
       #  units = "cm", 
       # width=10, height=10, 
       path="./output")

#### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- #### #### ----- ----- perform the comparisons
## import TPM table
TPM <- read.delim("./data/my_RNAseq_data/Ectocarpus-sp7/salmon.merged.gene_tpm.tsv", header = T)
colnames(TPM)[1] <- "ID"
# TPM <- TPM[,c(1,2,3,4,6,5,7,8)]
colnames(TPM)[c(3:6)] <- c("EcF1", "EcF2", "EcM1", "EcM2")
TPM <- TPM[,-2]
head(TPM)
# TPM$ID <- gsub('\\..*','',TPM$ID)
# TPM <- TPM %>%
#   group_by(ID) %>%  # Group by gene_id
#   summarise(across(everything(), sum, na.rm = TRUE)) # Sum all numeric columns

TPM$Female <- rowMeans(TPM[,2:3])
TPM$Male <- rowMeans(TPM[,4:5])
head(TPM)
#### summarise total genes detectable
x <- TPM
rownames(x) <- x$ID
x <- x[,-1]
row_sub = apply(x, 1, function(row) all(row != 0 ))
TPM.detected <- TPM[row_sub,]

#### perform the DESEQ comparisons
contrasts <- as.data.frame(results(run.dds, contrast = c("genotype","Male", "Female")))  ## use <lfcShrink> instead to give the MLE log2FC
contrasts$ID <- rownames(contrasts)
contrasts <- merge(contrasts, TPM, by = "ID")


#### VOLCANO PLOTS
#pdf("/ebio/abt5_projects/ceramiales/", width = 8, height = 10)

volcano.data <- contrasts
volcano.data$colour <- "dark gray"
volcano.data$colour[volcano.data$`log2FoldChange` < -1 & volcano.data$`padj` <= 0.05] <- "#D73027"
volcano.data$colour[volcano.data$`log2FoldChange` > 1 & volcano.data$`padj` <= 0.05] <- "#4575B4"
volcano.data$bias <- "Not significant"
volcano.data$bias[volcano.data$`log2FoldChange` < -1 & volcano.data$`padj` <= 0.05] <- "Female-biased genes"
volcano.data$bias[volcano.data$`log2FoldChange` > 1 & volcano.data$`padj` <= 0.05] <- "Male-biased genes"

# Create a data frame with the subset information
subset_data <- data.frame(
  direction = c("Female-biased genes", "Male-biased genes"),
  count = c(
    nrow(subset(volcano.data, log2FoldChange < -1 & padj <= 0.05)),
    nrow(subset(volcano.data, log2FoldChange > 1 & padj <= 0.05))
  )
)

# Create the ggplot
v <- ggplot(volcano.data, aes(x = log2FoldChange, y = -log10(padj), col = colour)) +
  geom_point(size = 0.5) +
  labs(
    title = "Differential gene expression analysis",
    x = "log2FoldChange",
    y = "-log10(adj pval)",
    caption = "Red genes are more expressed in the female and blue genes are more expressed in the male") +
  xlim(-20, 20) +
  ylim(0, 160) +
  geom_text(data = subset_data, aes(label = direction, x = c(-15, 15), y = 95), col = c("#D73027", "#4575B4"), size = 4) +
  geom_text(data = subset_data, aes(label = count, x = c(-15, 15), y = 85), col = c("#D73027", "#4575B4"), size = 5) +
  geom_hline(yintercept = 1.3, linetype = "dashed", colour = "#2a2727") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#2a2727") +
  scale_color_identity() +
  theme_minimal()

ggsave(filename = "volcanoplot_Ectocarpus-sp7-WT-lastassembly.pdf", 
       plot = last_plot(),
       #  units = "cm", 
       # width=10, height=10, 
       path="./output")

### ------------ plot TPM zscore of ALL DEGs

updown <- subset(contrasts, padj <= 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

matrix <- subset(TPM, ID %in% updown$ID)[,c(2:5)]
rownames(matrix) <- updown$ID

matrix <- t(scale(t(matrix), center = T, scale = T))
colours <- rev(brewer.pal(n = 9, name = "RdGy"))

pdf("./output/02_volcano_plots/plotheatmap_Ectocarpus-sp7-WT-lastassembly.pdf", width = 6, height = 6)
pheatmap(matrix, 
         cluster_rows = T,
         clustering_method = "ward.D2",
         cluster_cols = T,
         show_rownames = F,
         scale = "none",
         border_color = NA, 
         color = colours,
         #breaks = breaksList,
         fontsize = 12)
dev.off()

head(volcano.data)
print(v)

## Step 2: Ensure volcano.data is available -------------------------------------
if (!exists("volcano.data")) { 
  stop("The script did not generate volcano.data. Please check the script.")
}
## Step 3: Merge volcano.data with the master table by gene_id -------------------
# Rename the column gene_id by trad_id
names(volcano.data)[1] <- 'gene_id'
volcano.data <- volcano.data %>%  select(gene_id, bias) %>% 
  mutate(bias = ifelse(bias == 'Not significant', NA, bias)) %>% 
  distinct()
master_table <- left_join(master_table, volcano.data, by = "gene_id")
head(master_table %>%  select(gene_id, TPM.Female, TPM.Male, bias))
```

# ---------------------------------------------
# BLOCK5: Location on genome: sex chr / PAR / autosomes
# ---------------------------------------------

```{r}
master_table <- master_table %>%
  mutate(
    location = case_when(
      seq_id == 'Ec25_SDR_F' ~ 'Female-SDR',
      seq_id == 'chr_13' & gene_start >= 2235561 
      & gene_end <= 3158904 ~ 'Male-SDR',
      seq_id == 'chr_13' ~ 'PAR',
      TRUE ~ 'autosome'
    )
  )
head(master_table %>% select(gene_id, seq_id, location))
head(master_table %>% select(gene_id, seq_id, location) %>% filter(location == 'autosome'))
head(master_table %>% select(gene_id, seq_id, location) %>% filter(location == 'Male-SDR'))
head(master_table %>% select(gene_id, seq_id, location) %>% filter(location == 'Female-SDR'))
head(master_table %>% select(gene_id, seq_id, location) %>% filter(location == 'PAR'))
```

# ---------------------------------------------
# BLOCK6: Quartiles of expressed genes
# ---------------------------------------------

```{r, warning=FALSE}

# Define function to assign quantiles and write BED files
define_quantiles_and_write_bed <- function(df, tpm_column, sex, output_dir) {
  
  # Categorize expression
  df <- df %>%
    mutate(quantile = ifelse(.data[[tpm_column]] < 2, 'no_expression', 'expression'))
  
  # Compute summary statistics for 'expression' group
  summary_stats <- df %>%
    filter(quantile == 'expression') %>%
    dlookr::describe(.data[[tpm_column]], quantiles = c(0, 0.25, 0.5, 0.75))
  
  # Extract quartiles
  q0 <- summary_stats %>% pull(p00)
  q25 <- summary_stats %>% pull(p25)
  q50 <- summary_stats %>% pull(p50)
  q75 <- summary_stats %>% pull(p75)
  
  # Ensure quartiles are valid
  if (any(is.na(c(q25, q50, q75)))) {
    stop(paste("Error: Quartile values could not be extracted for", sex))
  }
  
  # Assign quantile groups
  df <- df %>%
    mutate(!!paste0("quantiles.", sex) := case_when(
      .data[[tpm_column]] < 2 ~ 'no_expression',
      .data[[tpm_column]] >= 2 & .data[[tpm_column]] <= q25 ~ '1st_quantile',
      .data[[tpm_column]] > q25 & .data[[tpm_column]] <= q50 ~ '2nd_quantile',
      .data[[tpm_column]] > q50 & .data[[tpm_column]] <= q75 ~ '3rd_quantile',
      .data[[tpm_column]] > q75 ~ '4th_quantile'
    ))
  df <- df %>% select(-quantile)
  # Update original dataframe
  assign("master_table", df, envir = .GlobalEnv)
  
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Function to write BED files
  write_bed <- function(df, quantile_name) {
    df_filtered <- df %>% filter(.data[[paste0("quantiles.", sex)]] == quantile_name)
    if (nrow(df_filtered) > 0) {
      df_filtered <- df_filtered %>% arrange(seq_id) %>% 
        select(seq_id, gene_start, gene_end)
      df_filtered$seq_id <- gsub("_", "", df_filtered$seq_id)  # Ensure consistency
      fwrite(df_filtered, file = paste0(output_dir, "/", quantile_name, 
                                        "_genes_", sex, ".Ectocarpus-sp7.bed"),
             sep = "\t", col.names = FALSE)
    }
  }
  
  # Write BED files for each quantile
  unique(df[[paste0("quantiles.", sex)]]) %>%
    purrr::walk(~ write_bed(df, .x))
}

# Set output directory
output_bed <- './data/05_OverlapEnrichment/Ectocarpus-sp7/'

# Apply function for Male and Female
define_quantiles_and_write_bed(master_table, "TPM.Female", "Female", output_bed)
define_quantiles_and_write_bed(master_table, "TPM.Male", "Male", output_bed)

print(master_table %>% select(gene_id, TPM.Female, TPM.Male, quantiles.Female, quantiles.Male))
```

# ---------------------------------------------
# BLOCK7: Chromatin states from hiHMM
# ---------------------------------------------

```{r, warning=FALSE}
# Load Female hiHMM file
hiHMM_female <- read.table("./hiHMM/hihmm.model2.K27.Ectocarpus-sp7-Ec25-FEMALE.recoloured.ReMapped.bed", header = FALSE, sep = '\t')

# Load Male hiHMM file
hiHMM_male <- read.table("./hiHMM/hihmm.model2.K27.Ectocarpus-sp7-Ec32-MALE.recoloured.ReMapped.bed", header = FALSE, sep = '\t')

# Modify seq_id as underscore is not tolerated by hiHMM
hiHMM_female$V1 <- gsub('chr', 'chr_', hiHMM_female$V1)
hiHMM_female$V1 <- gsub('HiCscaffold', 'HiC_scaffold_', hiHMM_female$V1)
hiHMM_female$V1 <- gsub('Ec25SDRF', 'Ec25_SDR_F', hiHMM_female$V1)
hiHMM_male$V1 <- gsub('chr', 'chr_', hiHMM_male$V1)
hiHMM_male$V1 <- gsub('HiCscaffold', 'HiC_scaffold_', hiHMM_male$V1)
hiHMM_male$V1 <- gsub('Ec25SDRF', 'Ec25_SDR_F', hiHMM_male$V1)

# Modify Emission labels
hiHMM_female$V4 <- paste0('E', hiHMM_female$V4)
hiHMM_male$V4 <- paste0('E', hiHMM_male$V4)

# Convert to GenomicRanges
hiHMM_gr_female <- GRanges(seqnames = hiHMM_female$V1,
                            ranges = IRanges(start = hiHMM_female$V2, 
                                             end = hiHMM_female$V3),
                            Emission = hiHMM_female$V4)

hiHMM_gr_male <- GRanges(seqnames = hiHMM_male$V1,
                          ranges = IRanges(start = hiHMM_male$V2, 
                                           end = hiHMM_male$V3),
                          Emission = hiHMM_male$V4)

mt_gr <- GRanges(seqnames = master_table$seq_id,
                          ranges = IRanges(start = master_table$gene_start, 
                                           end = master_table$gene_end),
                          geneid = master_table$gene_id)

# Find overlaps between genes and Emissions
overlaps_female <- findOverlaps(mt_gr, hiHMM_gr_female, minoverlap = 10)
overlaps_male <- findOverlaps(mt_gr, hiHMM_gr_male, minoverlap = 10)

# Extract emission states per gene for Female
emission_list_female <- tapply(hiHMM_gr_female$Emission[subjectHits(overlaps_female)], 
                               queryHits(overlaps_female), 
                               function(x) paste(unique(x), collapse = ","))

# Extract emission states per gene for Male
emission_list_male <- tapply(hiHMM_gr_male$Emission[subjectHits(overlaps_male)], 
                             queryHits(overlaps_male), 
                             function(x) paste(unique(x), collapse = ","))

# Initialize columns in master_table
master_table$Emissions.Female <- NA
master_table$Emissions.Male <- NA

# Assign values to corresponding rows
master_table$Emissions.Female[as.numeric(names(emission_list_female))] <- emission_list_female
master_table$Emissions.Male[as.numeric(names(emission_list_male))] <- emission_list_male

print(master_table %>% select(gene_id, Emissions.Female, Emissions.Male))
```

# ---------------------------------------------
# BLOCK8: Merge Emissions into chromatin Signature
# ---------------------------------------------

```{r}

melt_emissions_into_signature <- function(mt, sex) {
  # mt <- master_table 
  # mt <- master_table %>%  select(-c(Signature.Female, Signature.Male))
  # sex <- 'Female'
  mt <- mt %>% rename(Emissions = paste0('Emissions.', sex))
  patterns <- mt %>% select(gene_id, Emissions)
  states <- c(0:27)  # Define state range

  ## Summarize into a contingency table -------------------------------------
  ct <- patterns %>%
    mutate(Emissions = strsplit(as.character(Emissions), ",")) %>%
    unnest(Emissions)
  
  ct <- table(ct[, c("gene_id", "Emissions")])
  ct <- as.data.frame.matrix(ct)
  ct <- ct[,order(match(colnames(ct), states))]
  ct$gene_id <- rownames(ct)
  ct.filtered <- ct
  
  ## Handle Missing Columns ---------------------------------------------------
    expected_columns <- paste0('E',1:27)
    missing_cols <- setdiff(expected_columns, colnames(ct))
    
    if (length(missing_cols) > 0) {
      message(paste("Missing columns:", paste(missing_cols, collapse = ", "), "in Ectocarpus", sex))
      ct[missing_cols] <- 0
    }
    
    ## Group Emissions --------------------------------------------------------
  ct.filtered$TSS <- rowSums(ct[,c('E5','E10','E12')])
  ct.filtered$TSS <- ifelse(ct.filtered$TSS >= 1, 1, 0)  
  ct.filtered$GENE <- rowSums(ct[,c('E8','E2','E7','E9','E25')])
  ct.filtered$GENE <- ifelse(ct.filtered$GENE >= 1, 1, 0)  
  ct.filtered$M <- rowSums(ct[,c('E14','E17','E19','E4','E24','E18','E15','E11','E6','E20','E26','E22','E23')])
  ct.filtered$M <- ifelse(ct.filtered$M >= 1, 1, 0) 
  ct.filtered$pH <- rowSums(ct[,c('E27','E1')])
  ct.filtered$pH <- ifelse(ct.filtered$pH >= 1, 1, 0)
  ct.filtered$LS <- ifelse(rowSums(ct.filtered[,c('TSS','GENE','M','pH')] == 0) == 4, 1, 0)
  ct.filtered <- ct.filtered[,c("gene_id",'TSS','GENE','M','pH','LS')]
  ct.filtered$Signatures <- paste(ct.filtered[,2], ct.filtered[,3], ct.filtered[,4], ct.filtered[,5], ct.filtered[,6])
  length(unique(ct.filtered$Signatures))
  
  # Make States table ---------------------------------------------------
  states.matrix <- as.data.frame(names(table(ct.filtered$Signatures)))
  colnames(states.matrix) <- "combos"
  states.matrix <- separate(data = states.matrix, col = combos, sep = " ", into = c('TSS','GENE','M','pH','LS'), remove = T)
  rownames(states.matrix) <- paste("S", 1:nrow(states.matrix), sep = "")
  states.matrix$State <- paste("S", 1:nrow(states.matrix), sep = "")
  write.table(states.matrix, paste0("./output/merge_CS/hiHMM_signatures_05032025_Ectocarpus",sex,".txt"), row.names = T, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # Make Signatures table ---------------------------------------------------
  states.matrix$Signatures <- paste(states.matrix[,1], states.matrix[,2], states.matrix[,3], states.matrix[,4], states.matrix[,5]) 
  states.matrix
  head(ct.filtered)


  # Save the states table
  write.table(states.matrix, "./output/merge_CS/hiHMM_signatures_05032025_Ectocarpus_.txt", 
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

  # Merge with master table
  mt <- left_join(mt, select(ct.filtered, gene_id, Signatures), by = 'gene_id')
  mt <- left_join(mt, select(states.matrix, State, Signatures), by = 'Signatures')
  
  # Warning for missing state values
  missing_states <- sum(is.na(mt$State))
  if (missing_states > 0) {
    warning(paste("There are", missing_states, "NA values in State column. Check for missing emissions data."))
  }
  
  # Rename final column amd clean a bit
  mt <- mt %>% select(-Signatures)
  signature_col <- paste0('Signature.', sex)
  mt <- mt %>% rename(!!signature_col := State)  
  mt <- mt %>% rename(!!paste0('Emissions.', sex) := Emissions) 

  return(mt)
}
master_table <- melt_emissions_into_signature(master_table,'Female')
master_table <- melt_emissions_into_signature(master_table,'Male')
head(master_table %>%  select(gene_id, Emissions.Female, Emissions.Male, Signature.Female, Signature.Male))
head(master_table %>%  select(gene_id, Emissions.Female, Emissions.Male, Signature.Female, Signature.Male) %>% filter(Emissions.Female == Emissions.Male))
print(master_table %>%  select(gene_id, location, Emissions.Female, Emissions.Male, Signature.Female, Signature.Male) %>% filter(Emissions.Female == Emissions.Male) %>%  filter(Signature.Female != Signature.Male))
```

# ---------------------------------------------
# BLOCK9 - Bonus : Rename Signature depending on the order 
# ---------------------------------------------

```{r, eval = FALSE}
# # Import color signature table
# color_signature_table <- read_delim('./color_signature_table.csv', show_col_types = FALSE)
# 
# # Define signature columns to process
# signature_cols <- c("Signature.Female", "Signature.Male")
# 
# # Loop through each signature column
# for (signature_col in signature_cols) {
#   master_table <- master_table %>% rename(Signature = !!sym(signature_col))
#   master_table <- master_table %>%
#     left_join(color_signature_table %>% select(Signature, Signature_renamed), 
#               by = "Signature") %>% 
#     mutate(Signature := coalesce(Signature_renamed, Signature)) %>%  
#     select(-Signature_renamed) %>%   # Remove extra column after renaming 
#     rename(!!sym(signature_col) := Signature) # Rename to original
# }
# head(master_table %>%  select(gene_id, Emissions.Female, Emissions.Male, Signature.Female, Signature.Male))

```

# ---------------------------------------------
# BLOCK10 - Bonus : Plot Signature matrix 
# ---------------------------------------------

```{r}
# Import data
states.matrix <- read.delim("./output/merge_CS/hiHMM_signatures_05032025_Ectocarpus_Female.txt", header = TRUE, sep = "\t")

# Define order of states
x <- c(paste0("S", 16:1))

# # Rename and merge with color_signature_table
# states.matrix <- states.matrix %>%
#     rename(Signature = State) %>% 
#     left_join(color_signature_table %>% select(Signature, Signature_renamed), by = "Signature") %>% 
#     mutate(Signature = coalesce(Signature_renamed, Signature)) %>%  # Replace NAs with original Signature
#     select(-c(Signature_renamed, count)) %>%   
#     rename(State = Signature)  %>% # Rename back
#     arrange(match(State, x))

# Ensure only numeric columns are used for heatmap
numeric_matrix <- states.matrix %>% 
    select_if(is.numeric) %>%   # Select only numeric columns
    as.matrix()  # Convert to matrix

# Define the color scale
grey_scale <- colorRampPalette(c("white", "black"))(100)

# Open the PDF device
pdf("./output/merge_CS/SignaturesMatrix_05032025_Ectocarpus_Female.pdf", height = 20, width = 3)

# Create heatmap
s <- pheatmap(numeric_matrix,  
         cluster_rows = FALSE,        
         cluster_cols = FALSE,        
         show_rownames = TRUE,        
         show_colnames = TRUE,        
         scale = "none",              
         border_color = "white",      
         color = grey_scale, 
         fontsize = 10,               
         cellwidth = 20,              
         cellheight = 20)             

# Close PDF device
dev.off()


print(s)
```

### Final steps
```{r}
# Cleaning final master table format -----------------------------------------
master_table <- master_table %>% select(gene_id, gene_start, gene_end, seq_id, strand,
                      TPM.Female, log2_TPM_plus_1.Female, 
                      TPM.Male, log2_TPM_plus_1.Male, 
                      peak_K79.Female, peak_K79.Male,
                      quantiles.Female, quantiles.Male, 
                      bias, location,
                      Emissions.Female, Emissions.Male,
                      Signature.Female, Signature.Male)

# Correct grouping and NA assignment
master_table <- master_table %>% 
  mutate(
    TPM.Female = ifelse(location == 'Male-SDR', NA, TPM.Female),
    log2_TPM_plus_1.Female = ifelse(location == 'Male-SDR', NA, log2_TPM_plus_1.Female),
    Emissions.Female = ifelse(location == 'Male-SDR', NA, Emissions.Female),
    Signature.Female = ifelse(location == 'Male-SDR', NA, Signature.Female),
    peak_K79.Female = ifelse(location == 'Male-SDR', NA, peak_K79.Female),
    quantiles.Female = ifelse(location == 'Male-SDR', NA, quantiles.Female),
    Signature.Female = ifelse(location == 'Male-SDR', NA, Signature.Female),

    TPM.Male = ifelse(location == 'Female-SDR', NA, TPM.Male),
    log2_TPM_plus_1.Male = ifelse(location == 'Female-SDR', NA, log2_TPM_plus_1.Male),
    Emissions.Male = ifelse(location == 'Female-SDR', NA, Emissions.Male),
    Signature.Male = ifelse(location == 'Female-SDR', NA, Signature.Male),
    peak_K79.Male = ifelse(location == 'Female-SDR', NA, peak_K79.Male),
    quantiles.Male = ifelse(location == 'Female-SDR', NA, quantiles.Male),
    Signature.Male = ifelse(location == 'Female-SDR', NA, Signature.Male)
    )

print(master_table)

# Save the Final Master Table -----------------------------------------
write.table(master_table, file = './output/07_FINAL_MASTER_TABLES_w-Signature/fmts_Ectocarpus-sp7-longest-mRNA.csv', 
            sep = ';', col.names = TRUE, quote = FALSE, row.names = FALSE)

```
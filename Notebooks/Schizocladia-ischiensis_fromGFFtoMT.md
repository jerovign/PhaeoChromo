# MASTER table for Schizocladia ischiensis

### Preface
This notebook aims to keep track of how I made the MASTER TABLE that is used for all downstream analyses.
*gtf* used: `Schizocladia-ischiensis.gtf`, converted by AGAT as the original file is `Schizocladia-ischiensis.gff`. I had issue with the rtracklayer parser here, I am using the gff version. I had to use the gtf version for nf-core.

# Master table 

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
library(DESeq2)  # For variance-stabilizing transformation (VST)
library(ggplot2)

```

# ---------------------------------------------
# BLOCK1: Fill in gene_id gene_start gene_end with gff file 
# ---------------------------------------------

```{r}
## Extract gene_id, gene_start, and gene_end from gff file -------------------

# Define file path
gff_file <- "./data/gff_files/Schizocladia-ischiensis.gtf"
# Import GFF data
gff_data <- readGFF(gff_file)

# Filter for mRNA features
gff_data$type <- as.character(gff_data$type)
mRNA <- gff_data[gff_data$type == "transcript", ]

# Convert to GRanges object
mRNA_gr <- GRanges(
    seqnames = mRNA$seqid,
    ranges = IRanges(start = mRNA$start, end = mRNA$end),
    strand = mRNA$strand,  
    gene_id = mRNA$gene_id,
    trad_id = mRNA$transcript_id
)

## Fill in the values in the master table -----------------------------------
master_table <- as.data.table(mRNA_gr)
colnames(master_table) <- c("seq_id", "gene_start", "gene_end", 
                            "width", "strand", "gene_id", "trad_id")
print(master_table)
```


# ---------------------------------------------
# BLOCK2: Fetch TPM from salmon data
# ---------------------------------------------

# Quick check of the samples it seems there are 2 x triplicates.
```{r, warning=FALSE}
## --- read in counts and info
countsToUse <- read.delim("./data/my_RNAseq_data/Schizocladia-ischiensis/salmon.merged.gene_counts.tsv", header = T, as.is = T)
colnames(countsToUse)[c(3:8)] <- c("Si1", "Si2", "Si3", "Si4", "Si5", "Si6")
rownames(countsToUse) <-countsToUse$gene_id
countsToUse <- countsToUse[,-c(1,2)]
countsToUse <- as.matrix(countsToUse)
countsToUse <- subset(countsToUse, rowMax(countsToUse) >= 10)   #### this is optional

colData <- data.frame(triplicate = c(rep("A",3), rep("B",3)), 
                      row.names = colnames(countsToUse) )

dds <- DESeqDataSetFromMatrix(round(countsToUse), colData = colData, design=~triplicate)
dds <- estimateSizeFactors(dds)
run.dds <- DESeq(dds)

rld <- rlog(run.dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("triplicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$genotype <- factor(pcaData$triplicate, c("A","B"))
plotPCA <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=3, aes(color=triplicate)) +
  scale_color_manual(values =  c("#1f78b4", "#fb8072")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave(filename = "./data/my_RNAseq_data/Schizocladia-ischiensis/plotPCA_RNAseq_Schizocladia-ischiensis.pdf", 
       plot = plotPCA)
print(plotPCA)
```

```{r}
## Step 1: Read the Salmon output --------------------------------------------
salmon_file <- "./data/my_RNAseq_data/Schizocladia-ischiensis/salmon.merged.gene_tpm.tsv" 
salmon_data <- read.table(salmon_file, header = TRUE, sep = "\t")
colnames(salmon_data)[1] <- 'gene_id'

## Step 2: Calculate the mean to create the TPM column ----
salmon_data$TPM <- rowMeans(salmon_data[, 3:8], na.rm = TRUE)

## Step 3: Calculate log2(TPM + 1) and fill in the respective column ---------------
salmon_data$log2_TPM_plus_1 <- log2(salmon_data$TPM + 1)

salmon_data <- salmon_data %>% select(gene_id, TPM, log2_TPM_plus_1)
head(salmon_data)

# ## Step 4: Translate trad_id to gene_id thanks to a dictionnary (made with Python script not here) ----------------------- 

## Step 5: Merge TPM values with the master table based on trad_id -----------
master_table <- left_join(master_table, salmon_data, by = "gene_id")
head(master_table)
## Verify the merge was successful and there are no NA values where there shouldn't be ------
if (any(is.na(master_table$TPM))) {
  warning("Some geneIDs in the master table do not have matching TPM values in the Salmon output.")
}
```

# ---------------------------------------------
# BLOCK3: K79 peak annotation
# ---------------------------------------------

```{r, warning=FALSE}
# Load ChIP-seq peak file
peaks <- read.table("./data/peak_K79/Schizocladia-ischiensis-KU333_H3K79me2_merged-peaks_annotatePeaks.txt", header = TRUE, sep = '\t')
peaks_gr <- GRanges(seqnames = peaks$Chr,
                    ranges = IRanges(start = peaks$Start, 
                                     end = peaks$End
                                     ))

# Find overlaps between genes and peaks
overlaps <- findOverlaps(mRNA_gr, peaks_gr)

# Create a new column in the master table to indicate if a peak overlaps
master_table$peak_K79 <- ifelse(seq_along(master_table$gene_id) %in% queryHits(overlaps), TRUE, FALSE)

head(master_table)
```

# ---------------------------------------------
# BLOCK4: Sex-bias
# ---------------------------------------------

```{r}
master_table$bias <- 'NA'
```

# ---------------------------------------------
# BLOCK5: Location on genome: sex chr / PAR / autosomes
# ---------------------------------------------

```{r}
master_table$location <- 'NA'
```

# ---------------------------------------------
# BLOCK6: Quartiles of expressed genes
# ---------------------------------------------

```{r}
# Add 'quantiles' column based on TPM values
master_table <- master_table %>%
  mutate(quantiles = ifelse(TPM < 0.2, 'no_expression', 'expression'))

# Group by 'quantiles' and describe TPM
summary_stats <- master_table %>%
  group_by(quantiles) %>%
  dlookr::describe(TPM, quantiles = c(0, 0.25, 0.5, 0.75))

# Extract quartile values (make sure column names are correct)
q0 <- summary_stats %>% filter(quantiles == "expression") %>% pull(p00)
q25 <- summary_stats %>% filter(quantiles == "expression") %>% pull(p25)
q50 <- summary_stats %>% filter(quantiles == "expression") %>% pull(p50)
q75 <- summary_stats %>% filter(quantiles == "expression") %>% pull(p75)

# Ensure quartiles are not NULL
if (any(is.na(c(q25, q50, q75)))) {
  stop("Error: Quartile values could not be extracted. Check `summary_stats` output.")
}

# Reassign quantiles
master_table <- master_table %>%
  mutate(quantiles = case_when(
    TPM < 0.2 ~ 'no_expression',
    TPM >= 0.2 & TPM <= q25 ~ '1st_quantiles',
    TPM > q25 & TPM <= q50 ~ '2nd_quantiles',
    TPM > q50 & TPM <= q75 ~ '3rd_quantiles',
    TPM > q75 ~ '4th_quantiles'
  ))

# Ensure output directory exists
# output_bed <- './data/05_OverlapEnrichment/Schizocladia-ischiensis/'
output_bed <- './data/bed_files_for_plotProfile_quantiles/Schizocladia-ischiensis_KU333/'
dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)

# Function to write BED files
write_bed <- function(df, quantiles_name) {
  df <- df %>%
    arrange(seq_id) %>%
    select(seq_id, gene_start, gene_end, strand) ### need strand ?
  
  # df$seq_id <- gsub("_", "", df$seq_id)  # Ensure consistency
  
  fwrite(df, file = paste0(output_bed, "/", quantiles_name, "_genes.Schizocladia-ischiensis.bed"),
         sep = "\t", col.names = FALSE)
}

# Write BED files for each quantiles
write_bed(master_table %>% filter(quantiles == 'no_expression'), "no_expression")
write_bed(master_table %>% filter(quantiles == '1st_quantiles'), "1st_quantiles")
write_bed(master_table %>% filter(quantiles == '2nd_quantiles'), "2nd_quantiles")
write_bed(master_table %>% filter(quantiles == '3rd_quantiles'), "3rd_quantiles")
write_bed(master_table %>% filter(quantiles == '4th_quantiles'), "4th_quantiles")

```

# ---------------------------------------------
# BLOCK7: Chromatin states from hiHMM
# ---------------------------------------------

```{r, warning=FALSE}
# Load hiHMM file
hiHMM <- read.table("./data/hihmm.model2.K27.Schizocladia-ischiensis-KU333.recoloured.ReMapped.bed", header = FALSE, sep = '\t')
# Modify seq_id as underscore are not tolerated by hiHMM
hiHMM$V1 <- gsub('contig','_contig',hiHMM$V1)
hiHMM$V4 <- paste0('E',hiHMM$V4)
hiHMM_gr <- GRanges(seqnames = hiHMM$V1,
                    ranges = IRanges(start = hiHMM$V2, 
                                     end = hiHMM$V3),
                    Emission = hiHMM$V4)

# Find overlaps between genes and Emissions
overlaps <- findOverlaps(mRNA_gr, hiHMM_gr, minoverlap = 10)

# Extract emission states per gene
emission_list <- tapply(hiHMM_gr$Emission[subjectHits(overlaps)], queryHits(overlaps), function(x) paste(unique(x), collapse = ","))

# Initialize column and assign values
master_table$Emissions <- NA
master_table$Emissions[as.numeric(names(emission_list))] <- emission_list

```

# ---------------------------------------------
# BLOCK8: Merge Emissions into chromatin Signature
# ---------------------------------------------

```{r}
patterns <- master_table %>% select(gene_id, Emissions)

# Summarise into a contingency table --------------------------------------
states <- c(0:27) # hihmm model2, K0 = 7, K = 27, with all samples 

## HERMAPRHRODITE ------------------------------------------------------------------
ct <- patterns %>% 
  mutate(Emissions = strsplit(as.character(Emissions), ",")) %>%
  unnest(Emissions)

ct <- table(ct[, c("gene_id", "Emissions")])
ct <- as.data.frame.matrix(ct)
ct <- ct[,order(match(colnames(ct), states))]
ct$gene_id <- rownames(ct)
head(ct)

ct.filtered <- ct
### Grouping Emissions that looks the same --------------------------------
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
write.table(states.matrix, "./output/merge_CS/hiHMM_signatures_05032025.txt", row.names = T, col.names = TRUE, quote = FALSE, sep = "\t")

# Make Signatures table ---------------------------------------------------
states.matrix$Signatures <- paste(states.matrix[,1], states.matrix[,2], states.matrix[,3], states.matrix[,4], states.matrix[,5]) 
states.matrix
head(ct.filtered)

# MAKING FINAL MASTER TABLE -----------------------------------------
master_table <- left_join(master_table, ct.filtered, by = 'gene_id')
master_table <- left_join(master_table, states.matrix, by = 'Signatures')
master_table <- master_table %>% dplyr::rename(Signature = State)
if (any(is.na(master_table$State))) {
  warning("NA values in State.")
}

# Cleaning final master table fmt -----------------------------------------
master_table <- master_table %>% select(gene_id, gene_start,gene_end,seq_id, strand,
                      TPM,log2_TPM_plus_1,peak_K79,
                      quantiles,
                      bias,location,
                      Emissions,
                      Signature)
head(master_table)

write.table(master_table, file = './output/07_FINAL_MASTER_TABLES_w-Signature/fmts_Schizocladia-ischiensis.csv', sep = ';', col.names = TRUE, quote = FALSE, row.names = FALSE)
```

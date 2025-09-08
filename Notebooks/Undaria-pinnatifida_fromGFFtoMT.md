# Master table of Undaria pinnatifida


This notebook aims to keep track of how I made the MASTER TABLE that is used for all downstream analyses.
*gff* used: `PUBLIC_Undaria-pinnatifida_MALE.jv3.gff` derived from `PUBLIC_Undaria-pinnatifida_MALE.gff` as there were some issue with nf-core/pipeline.

#### Keeping track of GFF modifications

```{bash, eval = FALSE}
13.02.2025 Jeromine Vigneau: I addded 'HiC_scaffold_' before the scaffold number + AGAT + gffread
awk 'BEGIN {OFS="\t"} {if ($1 ~ /^[0-9]+$/) $1="HiC_scaffold_"$1; print}' PUBLIC_Undaria-pinnatifida_MALE.gtf > PUBLIC_Undaria-pinnatifida_MALE.jv.gtf
agat_convert_sp_gxf2gxf.pl -g PUBLIC_Undaria-pinnatifida_MALE.jv.gtf -o PUBLIC_Undaria_MALE.jv2.gtf
gffread PUBLIC_Undaria-pinnatifida_MALE.jv2.gtf -T -o PUBLIC_Undaria-pinnatifida_MALE.jv3.gtf
sed -E 's/;(\s*)$/\1/' PUBLIC_Undaria-pinnatifida_MALE.jv3.gtf > PUBLIC_Undaria-pinnatifida_MALE.jv4.gtf
 awk '{
    match($0, /transcript_id "([^"]+)"/, tid)
    match($0, /gene_id "([^"]+)"/, gid)

    if ($3 == "gene" && gid[1] != "") {
        print $0 " ID=\"" gid[1] "\";"
    }
    else if (($3 == "transcript" || $3 == "exon") && tid[1] != "" && gid[1] != "") {
        print $0 " ID=\"" tid[1] "\"; Parent=\"" gid[1] "\";"
    }
    else {
        print $0
    }
}' PUBLIC_Undaria-pinnatifida_MALE.jv4.gtf > PUBLIC_Undaria-pinnatifida_MALE.jv5.gtf
 replace " ID" with "; ID"
 remove trailing space again sed -E 's/;(\s*)$/\1/' PUBLIC_Undaria-pinnatifida_MALE.jv5.gtf > PUBLIC_Undaria-pinnatifida_MALE.jv6.gtf
 rename jv6 >> jv3 to reuse the file in nf-core pipeline and resume
 agat_convert_sp_gxf2gxf.pl -g PUBLIC_Undaria-pinnatifida_MALE.jv3.gtf -o PUBLIC_Undaria-pinnatifida_MALE.jv3.gff
 agat_sp_add_introns.pl --gff /ebio/abt5_projects/PhaeoChromo/data/gtf-gff-bed-files/PUBLIC_Undaria-pinnatifida_MALE.jv3.gff --out /ebio/abt5_projects/PhaeoChromo/data/gtf-gff-bed-files/PUBLIC_Undaria-pinnatifida_MALE.jv7.gff
```


#### Set up
```{r setup}
# Setting up environment
library(here)  # Optional but recommended for better path management
knitr::opts_knit$set(root.dir = "/home/jeromine/Documents/Scripts_Rstudio/MASTER_TABLES")
knitr::opts_chunk$set(echo = TRUE)
```

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
```


# ---------------------------------------------
# BLOCK1: Fill in gene_id gene_start gene_end with gff file 
# ---------------------------------------------

```{r}
## Extract gene_id, gene_start, and gene_end from gff file -------------------

# Define file path
gff_file <- "./data/gff_files/PUBLIC_Undaria-pinnatifida_MALE.jv3.gff"
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
    gene_id = as.character(gsub("\"|'", "", mRNA$Parent)),
    trad_id = gsub('transcript_id ', '', gsub("\"|'", "", mRNA$ID))
)

## Fill in the values in the master table -----------------------------------
master_table <- as.data.table(mRNA_gr) %>%distinct()
colnames(master_table) <- c("seq_id", "gene_start", "gene_end", 
                            "width", "strand", "gene_id", "trad_id")

# Remove duplicates due to transcripts (I am taking the whole gene, but gene models were lacking in this annotation)
master_table <- master_table %>% select(-c(trad_id,width)) %>%
  group_by(gene_id) %>%
  mutate(gene_start = min(gene_start),
         gene_end = max(gene_end)) %>%
  ungroup() %>%  # Corrected: Remove grouping
  distinct(gene_id, .keep_all = TRUE)  # Corrected: Keep one row per gene_id

# Convert into a GRanges again as I need it for later overlaps
mRNA_gr <- GRanges(
    seqnames = master_table$seq_id,  # Use "seq_id" instead of "seqid"
    ranges = IRanges(start = master_table$gene_start, end = master_table$gene_end),  # Use "gene_start" and "gene_end"
    strand = master_table$strand, 
    gene_id = master_table$gene_id
)


print(master_table)

```


# ---------------------------------------------
# BLOCK2: Fetch TPM from salmon data
# ---------------------------------------------

```{r}
## Step 1: Read the Salmon output --------------------------------------------
salmon_file <- "./data/my_RNAseq_data/Undaria-pinnatifida/salmon.merged.gene_tpm.tsv" 
salmon_data <- read.table(salmon_file, header = TRUE, sep = "\t")
colnames(salmon_data)[1] <- 'gene_id'

## Step 2: Calculate the mean to create the TPM column ----
salmon_data$TPM.Female <- rowMeans(salmon_data[, 3:5], na.rm = TRUE)
salmon_data$TPM.Male <- rowMeans(salmon_data[, 6:8], na.rm = TRUE)

## Step 3: Calculate log2(TPM + 1) and fill in the respective column ---------------
salmon_data$log2_TPM_plus_1.Female <- log2(salmon_data$TPM.Female + 1)
salmon_data$log2_TPM_plus_1.Male <- log2(salmon_data$TPM.Male + 1)

salmon_data <- salmon_data %>% select(gene_id, TPM.Female, TPM.Male, log2_TPM_plus_1.Female, log2_TPM_plus_1.Male)
head(salmon_data)

# ## Step 4: Translate trad_id to gene_id thanks to a dictionnary (made with Python script not here) ----------------------- 
# dict_trad <- read_csv('./data/Undaria-pinnatifida-mRNA_dict2.csv', show_col_types = FALSE)
# colnames(dict_trad) <- c('chromid', 'start', 'end', 'gene_id','trad_id')
# head(dict_trad)

## Step 5: Merge TPM values with the master table based on trad_id -----------
master_table <- left_join(master_table, salmon_data, by = "gene_id")
head(master_table %>% select(seq_id, gene_id, TPM.Female, TPM.Male, log2_TPM_plus_1.Female, log2_TPM_plus_1.Male))
## Verify the merge was successful and there are no NA values where there shouldn't be ------
if (any(is.na(master_table$TPM.Female | master_table$TPM.Male))) {
  warning("Some geneIDs in the master table do not have matching TPM values in the Salmon output.")
}
```

# ---------------------------------------------
# BLOCK3: K79 peak annotation
# ---------------------------------------------
### FEMALE

```{r}
# Load ChIP-seq peak file
peaks <- read.table("./data/peak_K79/Undaria_pinnatifida_female_PHAEO216_H3K79me2_merged-peaks_annotatePeaks.txt", header = TRUE, sep = '\t')
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
```{r}
# Load ChIP-seq peak file
peaks <- read.table("./data/peak_K79/Undaria_pinnatifida_male_PHAEO218_H3K79me2_merged-peaks_annotatePeaks.txt", header = TRUE, sep = '\t')
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
## Step 1: Source the script that generates volcano.data ------------------------
## volcano.data is obtain from DESeq2, a script not included here.
source("/home/jeromine/Documents/Scripts_Rstudio/DESeq2/code/deseq_Undaria-pinnatifida.R")
load("/home/jeromine/Documents/Scripts_Rstudio/DESeq2/output/volcano.data_Undaria-pinnatifida.RData")
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
  mutate(bias = ifelse(bias == 'Not significant', NA, bias))
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
      seq_id == 'HiC_scaffold_23' & gene_start >= 13867553 
      & gene_end <= 27276646 ~ 'Male-SDR',
      seq_id == 'HiC_scaffold_23' ~ 'PAR',
      TRUE ~ 'autosome'
    )
  )
head(master_table %>% select(gene_id, seq_id, location))
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
        select(seq_id, gene_start, gene_end, strand) ## need strand ?
      # df_filtered$seq_id <- gsub("_", "", df_filtered$seq_id)  # Ensure consistency
      fwrite(df_filtered, file = paste0(output_dir, "/", quantile_name, 
                                        "_genes_", sex, ".Undaria-pinnatifida.bed"),
             sep = "\t", col.names = FALSE)
    }
  }
  
  # Write BED files for each quantile
  unique(df[[paste0("quantiles.", sex)]]) %>%
    purrr::walk(~ write_bed(df, .x))
}

# Set output directory
# output_bed <- './data/05_OverlapEnrichment/Undaria-pinnatifida/'
output_bed <- '/home/jeromine/Documents/Scripts_Rstudio/MASTER_TABLES/data/bed_files_for_plotProfile_quantiles/Undaria-pinnatifida'

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
hiHMM_female <- read.table("/home/jeromine/Documents/Scripts_Rstudio/hiHMM/hiHMM/phaeochromo_august2024/06_output_browns-and-friends_K07/Remapped/hihmm.model2.K27.Undaria-pinnatifida-female-PHAEO216.recoloured.ReMapped.bed", header = FALSE, sep = '\t')

# Load Male hiHMM file
hiHMM_male <- read.table("/home/jeromine/Documents/Scripts_Rstudio/hiHMM/hiHMM/phaeochromo_august2024/06_output_browns-and-friends_K07/Remapped/hihmm.model2.K27.Undaria-pinnatifida-male-PHAEO218.recoloured.ReMapped.bed", header = FALSE, sep = '\t')

# Modify seq_id as underscore is not tolerated by hiHMM
hiHMM_female$V1 <- gsub('HiCscaffold', 'HiC_scaffold_', hiHMM_female$V1)
hiHMM_male$V1 <- gsub('HiCscaffold', 'HiC_scaffold_', hiHMM_male$V1)

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
    ## Group Emissions --------------------------------------------------------
    ############ Version 05 March 2025
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


  # Save the states table
  write.table(states.matrix, "./output/merge_CS/hiHMM_signatures_05032025_Undaria_.txt", 
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
head(master_table %>%  select(gene_id, Emissions.Female, Emissions.Male, Signature.Female, Signature.Male) %>% filter(Emissions.Female == Emissions.Male) %>%  filter(Signature.Female != Signature.Male))
```

# ---------------------------------------------
# BLOCK9 - Bonus : Rename Signature depending on the order 
# ---------------------------------------------

```{r}
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
# # Import data
# states.matrix <- read.delim("./output/merge_CS/hiHMM_signatures_25022025_Undaria_Female.txt", header = TRUE, sep = "\t")
# 
# # Define order of states
# x <- c(paste0("S", 16:1))
# 
# # Rename and merge with color_signature_table
# states.matrix <- states.matrix %>%
#     rename(Signature = State) %>% 
#     left_join(color_signature_table %>% select(Signature, Signature_renamed), by = "Signature") %>% 
#     mutate(Signature = coalesce(Signature_renamed, Signature)) %>%  # Replace NAs with original Signature
#     select(-Signature_renamed) %>%   
#     rename(State = Signature)  %>% # Rename back
#     arrange(match(State, x))
# 
# # Ensure only numeric columns are used for heatmap
# numeric_matrix <- states.matrix %>% 
#     select_if(is.numeric) %>%   # Select only numeric columns
#     as.matrix()  # Convert to matrix
# 
# # Define the color scale
# grey_scale <- colorRampPalette(c("white", "black"))(100)
# 
# # Open the PDF device
# pdf("./output/merge_CS/SignaturesMatrix_25022025_Undaria_Female.pdf", height = 20, width = 3)
# 
# # Create heatmap
# s <- pheatmap(numeric_matrix,  
#          cluster_rows = FALSE,        
#          cluster_cols = FALSE,        
#          show_rownames = TRUE,        
#          show_colnames = TRUE,        
#          scale = "none",              
#          border_color = "white",      
#          color = grey_scale, 
#          fontsize = 10,               
#          cellwidth = 20,              
#          cellheight = 20)             
# 
# # Close PDF device
# dev.off()
# 
# 
# print(s)
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
write.table(master_table, file = './output/07_FINAL_MASTER_TABLES_w-Signature/fmts_Undaria-pinnatifida.csv', 
            sep = ';', col.names = TRUE, quote = FALSE, row.names = FALSE)

```
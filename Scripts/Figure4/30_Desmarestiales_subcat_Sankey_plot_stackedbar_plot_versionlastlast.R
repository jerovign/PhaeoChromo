####################################################################
## Sankey + stacked bars + chi stats in one go (reuses SameSubSignature)
####################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(tools)
  library(conflicted)
  library(showtext)
})
conflicted::conflicts_prefer(dplyr::filter)

# ---- Fonts --------------------------------------------------------
hl_path <- systemfonts::match_fonts("Helvetica-Light")$path
font_add("helv_light", hl_path)
showtext_auto()

# ---- Paths --------------------------------------------------------
input_dir   <- './output/07_FINAL_MASTER_TABLES_w-Signature/FMT_relabeled/with_age_phylorank/'
output_dir  <- './output/31_SubSignature_synteny_stats/'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Stacked bar PDF
out_pdf_female     <- file.path(output_dir, "stacked_100_barplot_SignaturesChanges_stats_DESM_female.pdf")
out_pdf_male     <- file.path(output_dir, "stacked_100_barplot_SignaturesChanges_stats_DESM_male.pdf")

# ---- Colors & levels ----------------------------------------------
color_table <- read.csv('./color_signature_table.csv')
color_table$Signature_relabeled <- factor(color_table$Signature_relabeled,
                                          levels = color_table$Signature_relabeled)
category_order <- c('active','mixed','repressive','null')
category_colors <- c(active = '#E43A30', mixed = '#9392C8', repressive = '#2A74BB', null = '#999999')

# For stacked bars
cat_order <- c("Female-biased genes", "Male-biased genes","aSC",
               "PAR", "autosome", "Total genes")

fills <- c(
  "TRUE"                = "#E9E9E9",  # light grey for 'kept same signature'
  "Female-biased genes" = "#E43A30",  # red
  "Male-biased genes"   = "#2A74BB",  # blue
  "aSC" = "#1B733B", 
  "PAR"                 = "#1B733B",  # green
  "autosome"            = "#FFB000",  # amber
  "Total genes"         = "#999999"   # dark grey
)

# ---- Helpers -------------------------------------------------------
safe_filename <- function(x) {
  x <- gsub('[^A-Za-z0-9_-]+', '_', x)
  gsub('_+$','',x)
}

p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)  ~ " ",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# ---- Main process loop (Sankey + counts for bars) ---------------

## -- Load OG table --
OG_dir <- './data/Orthofinder/Results_Dec13/Orthogroups/'
og_data <-
  readr::read_tsv(paste0(OG_dir, 'Orthogroups.tsv')) |>
  tidyr::drop_na() |>
  tidyr::pivot_longer(
    cols = -Orthogroup,
    names_to = "species",
    values_to = "gene_id"
  ) |>
  dplyr::filter(!is.na(gene_id)) |>
  tidyr::separate_rows(gene_id, sep = ",\\s*") |>
  dplyr::filter(gene_id != "") |>
  # filter isoforms in Ddud
  dplyr::filter(!stringr::str_detect(gene_id, "\\.t[0-9]+$") | stringr::str_detect(gene_id, "\\.t1$"))

# get single copy orthogroup for DEsmarestiales
  sp_pair <- c("Desmarestia-dudresnayi_BRAKER_proteins",
               "Desmarestia-herbacea_MALE_proteins")

og_sc.pre <- og_data %>%
  dplyr::filter(species %in% sp_pair) %>%
  dplyr::count(Orthogroup, species, name = "n") %>%
  dplyr::filter(n == 1) %>%                       # single-copy in each of the pair
  dplyr::count(Orthogroup, name = "n_sp") %>%
  dplyr::filter(n_sp == length(sp_pair)) %>%      # both species present
  dplyr::select(Orthogroup)

og_sc <- og_data %>%
  dplyr::filter(Orthogroup %in% og_sc.pre$Orthogroup,
                species %in% sp_pair)

OG <- og_sc %>%
  select(Orthogroup, species, gene_id) %>%
  distinct() %>%
  pivot_wider(names_from = species, values_from = gene_id)
rownames(OG) <- OG[[1]]
# Keep only the columns of interest and normalize geneIDs
OG$`Desmarestia-dudresnayi_BRAKER_proteins` <- sapply(strsplit(as.character(OG$`Desmarestia-dudresnayi_BRAKER_proteins`), '\\.'), `[`, 1)
OG$`Desmarestia-herbacea_MALE_proteins`   <- gsub('prot', 'mRNA', OG$`Desmarestia-herbacea_MALE_proteins`)

## -- Initialisation --
files <- list.files(path = input_dir, full.names = TRUE, pattern = 'Desmarestia')
stopifnot(length(files) > 0)

### Desmarestia dudresnayi
data_dud <- read.table(file = files[1], sep = ';', header = TRUE, check.names = FALSE, quote = "\"")
data_dud <- data_dud %>%
  filter(!is.na(gene_id), !is.na(Signature)) %>% 
  mutate(Signature = factor(Signature, levels = levels(color_table$Signature_relabeled)))
map_dd <- data_dud %>% dplyr::select(gene_id, Signature, location) %>% 
  rename(location_Ddud = location)

### Desmarestia herbacea
data_herb <- read.table(file = files[2], sep = ';', header = TRUE, check.names = FALSE, quote = "\"")

bars_all <- tibble() # to collect comparisons CategoryType/Category counts & proportions & p-values

# loop to handle both sexes
sexes <- c('Signature.Female','Signature.Male') 
for (sex_col in sexes) {
  sex <- gsub("Signature\\.", "", sex_col)
  comp_name <- paste0("Desmarestiales_Ddud_vs_", sex)
  data_herb[[sex_col]] <- factor(data_herb[[sex_col]], levels = levels(color_table$Signature_relabeled))

  # per-iteration filter (use the variable column)
  data_herb_i <- data_herb %>%
    dplyr::filter(!is.na(location), !is.na(.data[[sex_col]]))
  
  # build map from the filtered copy
  map_dh <- data_herb_i %>%
    dplyr::select(gene_id, all_of(sex_col), bias, location) %>%
    dplyr::rename(Signature_herb = all_of(sex_col),
                  bias_Dherb = bias,
                  location_Dherb = location) %>% 
    dplyr::filter(!is.na(Signature_herb))
  # if analyzing Female, drop herb rows in Male SDR
  map_dh <- map_dh %>% 
    { if (identical(sex_col, "Signature.Female") &&
          "location_Dherb" %in% names(.)) {
      dplyr::filter(., is.na(location_Dherb) | !location_Dherb %in% "Male-SDR")
    } else . } 
  
  df_OG_signature <- OG %>%
    # attach Ddud signature
    left_join(map_dd, by = c(`Desmarestia-dudresnayi_BRAKER_proteins` = "gene_id")) %>%
    rename(Ddud = Signature) %>%
    # attach Dherb signature
    left_join(map_dh, by = c(`Desmarestia-herbacea_MALE_proteins` = "gene_id")) %>%
    rename(Dherb = Signature_herb) %>%
    # build sub-signatures
    mutate(
      SubSignature.Ddud = case_when(
        Ddud %in% c("S16","S15","S14","S13","S12") ~ "active",
        Ddud %in% c("S11","S10","S9","S8","S7","S6","S5") ~ "mixed",
        Ddud %in% c("S4","S3","S2") ~ "repressive",
        TRUE ~ "null"
      ),
      SubSignature.Dherb = case_when(
        Dherb %in% c("S16","S15","S14","S13","S12") ~ "active",
        Dherb %in% c("S11","S10","S9","S8","S7","S6","S5") ~ "mixed",
        Dherb %in% c("S4","S3","S2") ~ "repressive",
        TRUE ~ "null"
      )
    ) %>%
    mutate(
      SubSignature.Ddud  = factor(SubSignature.Ddud,  levels = category_order),
      SubSignature.Dherb = factor(SubSignature.Dherb, levels = category_order),
      SameSubSignature   = SubSignature.Ddud == SubSignature.Dherb
    ) %>% 
    filter(!is.na(Dherb),!is.na(Ddud),!is.na(location_Dherb),)

  # ---------- Sankey ----------
  sankey_data <- df_OG_signature %>%
    count(SubSignature.Ddud, SubSignature.Dherb, SameSubSignature, name = "Freq")
  
  same_genes <- sankey_data %>%
    filter(SameSubSignature) %>%
    summarise(n = sum(Freq), .groups = "drop") %>%
    pull(n)
  same_genes <- ifelse(length(same_genes) == 0, 0, same_genes)
  
  total_genes <- sum(sankey_data$Freq)
  if (total_genes == 0) {
    message("No genes after filtering for ", comp_name, " – skipping plots.")
    next
  }
  same_subsignature_percentage <- 100 * same_genes / total_genes
  
  # Export flux tables
  flux_counts <- sankey_data %>%
    group_by(SubSignature.Ddud, SubSignature.Dherb) %>%
    summarise(Gene_Count = sum(Freq), .groups = "drop")
  write.csv(
    flux_counts,
    file.path(output_dir, paste0("Flux_table_", safe_filename(comp_name), ".csv")),
    row.names = FALSE
  )
  
  flux_matrix <- tidyr::pivot_wider(
    flux_counts,
    names_from  = SubSignature.Dherb,
    values_from = Gene_Count,
    values_fill = 0
  ) %>%
    arrange(match(SubSignature.Ddud, category_order))
  write.csv(
    flux_matrix,
    file.path(output_dir, paste0("Flux_matrix_", safe_filename(comp_name), ".csv")),
    row.names = FALSE
  )
  
  # Dynamic axis labels
  left_axis_label  <- "Ddud"
  right_axis_label <- paste0("D. herbacea (", sex, ")")
  
  sankey_plot <- ggplot(
    sankey_data,
    aes(axis1 = SubSignature.Ddud, axis2 = SubSignature.Dherb, y = Freq)
  ) +
    geom_alluvium(
      aes(fill = ifelse(SameSubSignature, "#EBEBEB", as.character(SubSignature.Ddud))),
      width = 1/12, alpha = 0.9
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 1/6, color = "white", alpha = 1, 
    ) +
    # geom_text(
    #   stat = "stratum", aes(label = after_stat(stratum)),
    #   family = "helv_light", size = 4
    # ) +
    # geom_text(
    #   aes(label = Freq), stat = "flow",
    #   size = 4, color = "black", hjust = -1,
    #   family = "helv_light"
    # ) +
    scale_x_discrete(limits = c(left_axis_label, right_axis_label), expand = c(0.1, 0.1)) +
    scale_fill_manual(values = c("#EBEBEB" = "#EBEBEB", category_colors)) +
    labs(
      title = paste0(
        comp_name, "\n",
        same_genes, " genes keeping signature ",
        sprintf("%.1f%%", same_subsignature_percentage)
      ),
      x = NULL, y = "Number of genes"
    ) +
    theme_void()+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14),
      axis.title  = element_blank(),
      axis.text.x = element_text(size = 14)
    )
  
  ggsave(
    file.path(output_dir, paste0("Synteny_Signature_4subcat_", safe_filename(comp_name), ".pdf")),
    plot = sankey_plot, width = 4, height = 5
  )
  
  # ---------- Build bars data from SameSubSignature ----------
  # Clean bias (drop "Not significant") if present
  data_for_bars <- df_OG_signature %>%
    mutate(bias_Dherb = if ("bias_Dherb" %in% names(.)) ifelse(bias_Dherb == "Not significant", NA, bias_Dherb) else NA)
  
  # Long format over available CategoryTypes
  long_cat <- data_for_bars %>%
    pivot_longer(
      cols = any_of(c("bias_Dherb", "location_Ddud")),
      names_to = "CategoryType", values_to = "Category"
    ) %>%
    filter(!is.na(Category))
  
  # Per CategoryType/Category × TRUE/FALSE counts
  ct_counts <- long_cat %>%
    group_by(CategoryType, Category, SameSubSignature) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Overall bar (single category)
  overall_counts <- data_for_bars %>%
    transmute(
      CategoryType = "Overall",
      Category     = "Total genes",
      SameSubSignature = SameSubSignature
    ) %>%
    group_by(CategoryType, Category, SameSubSignature) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Combine and compute proportions
  counts_all <- bind_rows(ct_counts, overall_counts) %>%
    mutate(SameSubSignature = ifelse(SameSubSignature, "TRUE", "FALSE")) %>%
    group_by(CategoryType, Category) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  ## 2×2 chi-square per bar vs Overall (Total genes) ----------------
  ov <- counts_all %>%
    filter(CategoryType == "Overall", Category == "Total genes") %>%
    summarise(
      ov_true  = sum(Count[SameSubSignature == "TRUE"],  na.rm = TRUE),
      ov_false = sum(Count[SameSubSignature == "FALSE"], na.rm = TRUE)
    ) %>% as.list()
  
  per_bar <- counts_all %>%
    filter(!(CategoryType == "Overall" & Category == "Total genes")) %>%
    group_by(CategoryType, Category) %>%
    summarise(
      n_true  = sum(Count[SameSubSignature == "TRUE"],  na.rm = TRUE),
      n_false = sum(Count[SameSubSignature == "FALSE"], na.rm = TRUE),
      .groups = "drop"
    )
  
  pvals <- per_bar %>%
    rowwise() %>%
    mutate(
      test            = list(suppressWarnings(
        chisq.test(rbind(c(n_true, n_false), c(ov$ov_true, ov$ov_false)), rescale.p = TRUE, simulate.p.value = TRUE, B = 1e6)
      )),
      chi_sq          = test$statistic,
      chi_df          = test$parameter,
      chi_raw_p_value = test$p.value
    ) %>%
    ungroup() %>%
    mutate(chi_adjusted_p_value = p.adjust(chi_raw_p_value, method = "bonferroni")) %>%
    select(CategoryType, Category, chi_raw_p_value, chi_adjusted_p_value)
  
  bars_species <- counts_all %>%
    left_join(pvals, by = c("CategoryType", "Category")) %>%
    mutate(Comparison = comp_name)
  
  # Save per-comparison CSV
  write.csv(
    bars_species,
    file.path(output_dir, paste0("table_barplot_w-chi-test_", safe_filename(comp_name), ".csv")),
    row.names = FALSE
  )
  
  # --- Make the stacked bar for THIS comparison (inside the loop) ---
  df_sp <- bars_species %>%
    mutate(
      Category = ifelse(Category == "Overall", "Total genes", Category),
      FillKey  = ifelse(SameSubSignature == "TRUE", "TRUE", as.character(Category)),
      Category = factor(Category, levels = cat_order)
    )
  
  # stars for this comparison
  df_sig <- df_sp %>%
    group_by(Category) %>%
    summarise(p = suppressWarnings(min(chi_adjusted_p_value, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(stars = p_to_stars(p))
  
  p <- ggplot(df_sp, aes(x = Category, y = Proportion, fill = FillKey)) +
    geom_col(width = 0.5, position = "stack") +
    geom_text(aes(label = Count),
              position = position_stack(vjust = 0.5),
              size = 4.2, family = "helv_light", angle = 90) +
    geom_text(
      data = df_sig,
      aes(x = Category, y = 1.02, label = stars),
      inherit.aes = FALSE,
      hjust = 0.4, vjust = -0.1,
      size = 4.5, family = "helv_light", fontface = "bold"
    ) +
    scale_fill_manual(values = fills, breaks = c("TRUE", cat_order), drop = FALSE) +
    guides(fill = "none") +
    theme_void() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 14, angle = 90, family = "helv_light"),
      axis.text.y = element_blank(),
      plot.margin = margin(t = 8, r = 30, b = 12, l = 12)
    )
  
  ggsave(file.path(output_dir, paste0("stacked_barplot_", safe_filename(comp_name), ".pdf")),
         plot = p, width = 5, height = 6, dpi = 300)
  message("Saved: stacked_barplot_", comp_name)
  
  # accumulate into one DF for the final combined CSV
  bars_all <- dplyr::bind_rows(bars_all, bars_species)

}
  
write.csv(bars_all,
          file.path(output_dir, "table_barplot_w-chi-test_Desmarestiales.csv"),
          row.names = FALSE)
message("Saved stacked chi-square table: table_barplot_w-chi-test_Desmarestiales.csv")

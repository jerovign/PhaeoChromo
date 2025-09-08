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
out_pdf     <- file.path(output_dir, "stacked_100_barplot_SignaturesChanges_stats.pdf")

# ---- Colors & levels ----------------------------------------------
color_table <- read.csv('./color_signature_table.csv')
color_table$Signature_relabeled <- factor(color_table$Signature_relabeled,
                                          levels = color_table$Signature_relabeled)

category_order <- c('active','mixed','repressive','null')
category_colors <- c(active = '#E43A30', mixed = '#9392C8', repressive = '#2A74BB', null = '#999999')

# For stacked bars
cat_order <- c("Female-biased genes", "Male-biased genes","aSC",
               "PAR", "autosome", "Total genes")

species_order <- c(
  "Ectocarpus-sp7-longest-mRNA",
  "Scytosiphon-promiscuus",
  "Undaria-pinnatifida",
  "Desmarestia-herbacea",
  "Desmarestia-dudresnayi"
)
species_labels <- c("Esp7","Spro","Upin","Dher","Ddud")
species_map    <- setNames(species_labels, species_order)

fills <- c(
  "TRUE"                = "#E9E9E9",  # light grey for 'kept same signature'
  "Female-biased genes" = "#E43A30",  # pink
  "Male-biased genes"   = "#2A74BB",  # blue
  "aSC" = 'violet', 
  "PAR"                 = "#1B733B",  # sand
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

# ---- Main processing loop (Sankey + counts for bars) ---------------
files <- list.files(path = input_dir, full.names = TRUE)
stopifnot(length(files) > 0)

bars_all <- list()  # to collect per-species CategoryType/Category counts & proportions & p-values

for (file in files) {
  species_name <- tools::file_path_sans_ext(basename(file)) |> gsub('fmts_', '', x = _)
  
  data <- tryCatch({
    read.table(file = file, sep = ';', header = TRUE, check.names = FALSE, quote = "\"")
  }, error = function(e) NULL)
  
  if (is.null(data)) {
    message('Skipping ', species_name, ' - cannot read file')
    next
  }
  if (!all(c('Signature.Female','Signature.Male') %in% colnames(data))) {
    message('Skipping ', species_name, ' - missing required columns')
    next
  }
  
  data <- data %>%
    filter(!is.na(gene_id), !is.na(Signature.Female), !is.na(Signature.Male))
  
  data$Signature.Female <- factor(data$Signature.Female, levels = levels(color_table$Signature_relabeled))
  data$Signature.Male   <- factor(data$Signature.Male,   levels = levels(color_table$Signature_relabeled))
  
  # Build SubSignatures and SameSubSignature (this will be reused for bars)
  data_to_plot <- data %>%
    select(gene_id, Signature.Female, Signature.Male, dplyr::any_of(c("bias","location"))) %>%
    mutate(
      SubSignature.Female = case_when(
        Signature.Female %in% c('S16','S15','S14','S13','S12') ~ 'active',
        Signature.Female %in% c('S11','S10','S9','S8','S7','S6','S5') ~ 'mixed',
        Signature.Female %in% c('S4','S3','S2') ~ 'repressive',
        TRUE ~ 'null'
      ),
      SubSignature.Male = case_when(
        Signature.Male %in% c('S16','S15','S14','S13','S12') ~ 'active',
        Signature.Male %in% c('S11','S10','S9','S8','S7','S6','S5') ~ 'mixed',
        Signature.Male %in% c('S4','S3','S2') ~ 'repressive',
        TRUE ~ 'null'
      )
    ) %>%
    mutate(
      SubSignature.Female = factor(SubSignature.Female, levels = category_order),
      SubSignature.Male   = factor(SubSignature.Male,   levels = category_order),
      SameSubSignature    = SubSignature.Female == SubSignature.Male
    )
  
  # ---------- Sankey ----------
  sankey_data <- data_to_plot %>%
    count(SubSignature.Female, SubSignature.Male, SameSubSignature, name = 'Freq')
  
  same_genes <- sankey_data %>% filter(SameSubSignature) %>% summarise(n = sum(Freq), .groups="drop") %>% pull(n)
  same_genes <- ifelse(length(same_genes) == 0, 0, same_genes)
  total_genes <- sum(sankey_data$Freq)
  if (total_genes == 0) {
    message('No genes after filtering for ', species_name, ' – skipping plots.')
    next
  }
  same_subsignature_percentage <- 100 * same_genes / total_genes
  
  # Export flux tables
  flux_counts <- sankey_data %>%
    group_by(SubSignature.Female, SubSignature.Male) %>%
    summarise(Gene_Count = sum(Freq), .groups = 'drop')
  write.csv(flux_counts,
            file.path(output_dir, paste0('Flux_table_', safe_filename(species_name), '.csv')),
            row.names = FALSE)
  
  flux_matrix <- tidyr::pivot_wider(flux_counts,
                                    names_from = SubSignature.Male,
                                    values_from = Gene_Count,
                                    values_fill = 0) %>%
    arrange(match(SubSignature.Female, category_order))
  write.csv(flux_matrix,
            file.path(output_dir, paste0('Flux_matrix_', safe_filename(species_name), '.csv')),
            row.names = FALSE)
  
  sankey_plot <- ggplot(sankey_data, aes(axis1 = SubSignature.Female, axis2 = SubSignature.Male, y = Freq)) +
    geom_alluvium(
      aes(fill = ifelse(SameSubSignature, "#EBEBEB", as.character(SubSignature.Female))),
      width = 1/12, alpha = 0.9
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 1/6, color = "white", alpha = 1, 
    ) +
    # geom_text(
    #   stat = 'stratum', aes(label = after_stat(stratum)),
    #   family = "helv_light", size = 4
    # ) +
    # geom_text(
    #   aes(label = Freq), stat = 'flow',
    #   size = 4, color = 'black', hjust = -1,
    #   family = "helv_light"
    # ) +
    scale_x_discrete(limits = c('Female', 'Male'), expand = c(0.1, 0.1)) +
    scale_fill_manual(values = c("#EBEBEB" = "#EBEBEB", category_colors)) +
    labs(
      title = paste0(
        species_name, '\n',
        same_genes, ' genes keeping signature ',
        sprintf('%.1f%%', same_subsignature_percentage)
      ),
      x = NULL, y = 'Number of genes'
    ) +
    theme_void()+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14),
      axis.title  = element_blank(),
      axis.text.x = element_text(size = 14)
    )
  
  ggsave(file.path(output_dir, paste0('Synteny_Signature_4subcat_', safe_filename(species_name), '.pdf')),
         plot = sankey_plot, width = 4, height = 5)
  
  # ---------- Build bars data from SameSubSignature ----------
  # Clean bias (drop "Not significant") if present
  data_for_bars <- data_to_plot %>%
    mutate(bias = if ("bias" %in% names(.)) ifelse(bias == "Not significant", NA, bias) else NA)
  
  # Long format over available CategoryTypes
  long_cat <- data_for_bars %>%
    pivot_longer(cols = any_of(c("bias","location")),
                 names_to = "CategoryType", values_to = "Category") %>%
    filter(!is.na(Category))
  
  # Per CategoryType/Category × TRUE/FALSE counts
  ct_counts <- long_cat %>%
    group_by(CategoryType, Category, SameSubSignature) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Overall bar (single category)
  overall_counts <- data_for_bars %>%
    transmute(CategoryType = "Overall",
              Category     = "Total genes",
              SameSubSignature = SameSubSignature ) %>%
    group_by(CategoryType, Category, SameSubSignature) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Combine and compute proportions
  counts_all <- bind_rows(ct_counts, overall_counts) %>%
    mutate(SameSubSignature = ifelse(SameSubSignature, "TRUE", "FALSE")) %>%
    group_by(CategoryType, Category) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  ## 2×2 chi-square per bar vs Overall (Total genes) ----------------
  # overall TRUE/FALSE (from counts_all)
  ov <- counts_all %>%
    filter(CategoryType == "Overall", Category == "Total genes") %>%
    summarise(
      ov_true  = sum(Count[SameSubSignature == "TRUE"],  na.rm = TRUE),
      ov_false = sum(Count[SameSubSignature == "FALSE"], na.rm = TRUE)
    ) %>% as.list()
  
  # per (CategoryType, Category) TRUE/FALSE totals (exclude Overall row)
  per_bar <- counts_all %>%
    filter(!(CategoryType == "Overall" & Category == "Total genes")) %>%
    group_by(CategoryType, Category) %>%
    summarise(
      n_true  = sum(Count[SameSubSignature == "TRUE"],  na.rm = TRUE),
      n_false = sum(Count[SameSubSignature == "FALSE"], na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2×2 chi-square vs overall; keep raw and Bonferroni-adjusted p
  pvals <- per_bar %>%
    rowwise() %>%
    mutate(
      test            = list(suppressWarnings(
        chisq.test(rbind(c(n_true, n_false),
                         c(ov$ov_true, ov$ov_false)))
      )),
      chi_sq          = test$statistic,
      chi_df          = test$parameter,
      chi_raw_p_value = test$p.value
    ) %>%
    ungroup() %>%
    mutate(chi_adjusted_p_value = p.adjust(chi_raw_p_value, method = "bonferroni")) %>%
    select(CategoryType, Category, chi_raw_p_value, chi_adjusted_p_value)
  
  
  # Assemble per-species bars table
  bars_species <- counts_all %>%
    left_join(pvals, by = c("CategoryType","Category")) %>%
    mutate(Species = species_name)
  
  # Save per-species CSV 
  write.csv(bars_species,
            file.path(output_dir, paste0("table_barplot_w-chi-test_", safe_filename(species_name), ".csv")),
            row.names = FALSE)
  
  bars_all[[species_name]] <- bars_species

# Loop over species in bars_all and plot each stacked bar separately
for (sp in names(bars_all)) {
  df_sp <- bars_all[[sp]] %>%
    mutate(
      Category = ifelse(Category == "Overall", "Total genes", Category),
      FillKey  = ifelse(SameSubSignature == "TRUE", "TRUE", as.character(Category)),
      Category = factor(Category, levels = cat_order)
    )
  
  # Stars for this species
  df_sig <- df_sp %>%
    group_by(Category) %>%
    summarise(p = suppressWarnings(min(chi_adjusted_p_value, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(stars = p_to_stars(p))
  
  p <- ggplot(df_sp, aes(x = Category, y = Proportion, fill = FillKey)) +
    geom_col(width = 0.5, position = "stack") +
    geom_text(aes(label = Count),
              position = position_stack(vjust = 0.5),
              size = 4.2, family = "helv_light") +
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
  
  ggsave(file.path(output_dir, paste0("stacked_barplot_", safe_filename(sp), ".pdf")),
         plot = p, width = 5, height = 6, dpi = 300)
  message("Saved: stacked_barplot_", sp)
}}

# # ---- Stack all per-species chi-square tables into one -------------
# all_tables <- list.files(output_dir,
#                          pattern = "^table_barplot_w-chi-test_.*\\.csv$",
#                          full.names = TRUE) |>
#   purrr::map_dfr(~ read.csv(.x) |>
#                    mutate(Species = gsub("^table_barplot_w-chi-test_", "",
#                                          tools::file_path_sans_ext(basename(.x)))),
#                  .id = NULL)
# 
# # Sauvegarde dans un fichier unique
# write.csv(all_tables,
#           file.path(output_dir, "table_barplot_w-chi-test_ALL_species.csv"),
#           row.names = FALSE)
# 
# message("Saved stacked chi-square table: table_barplot_w-chi-test_ALL_species.csv")

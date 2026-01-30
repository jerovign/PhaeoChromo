# How to run the scripts

1. Set this directory as a project on RStudio

2. Run the `.Rmd` files in this directory.

# Files in this directory

* `RNA-seq_correlation_heatmap`: correlation between RNA-seq samples used in this paper.

* `K79_expression`: expression levels of genes marked by H3K79me2 (for whole gene or just TSS)

* `signature_dist_sex_chromosome`: differences in signature between sex chromosomes and autosomes.

* `signature_conservation`: comparison between species to detect the number of orthologues that share signatures - compared to permutation results.

* `signature_conservation_sex`: signatures are more conserved between species than permutations. The next question is how this looks if we analyse genes in the sex chromosomes separately from autosomes.

* `single_copy_analysis`: analysis of genes in the one-to-one (single copy) orthologue set.

* `signature_phylorank_tau`: the distribution of gene age (phylorank) and expression dyanamics (tau) for each signature. This is done for all species.

* `signature_orphan`: number of orphan genes in each signature. This is done for all species.

> [!NOTE]
> `|>` pipes are used here (introduced in `R 4.1.0`). Otherwise, `|>` should be replaced with `magrittr`'s `%>%`. Furthermore, `duckplyr v0.4.1` was used but the package has been updated since, resulting in some backcompatibility issues. This is easily resolved by changing `duckplyr::` to `dplyr::`, except for `duckplyr::as_duckdb_tibble()`.

Further scripts used in the revisions are saved as PDFs in `revision_pdf`.

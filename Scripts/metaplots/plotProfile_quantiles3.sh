#!/bin/bash
#$ -pe parallel 6
#$ -l h_vmem=16G
#$ -j y
#$ -cwd

# Base directories
BASE_DIR="/ebio/abt5_projects/PhaeoChromo/data"
PROJECT_DIR="$BASE_DIR/early-stage_evochromo"
OUTPUT_DIR="$BASE_DIR/plotProfile/plots_per_mark"
mkdir -p "$OUTPUT_DIR"  # Ensure output directory exists

# Define samples and species
SAMPLES=(
    # "Desmarestia-dudresnayi_SxS43"
    # "Desmarestia-herbacea_SxS46_FEMALE"
    # "Desmarestia-herbacea_SxS47_MALE"
    # "Schizocladia-ischiensis_KU333"
    # "Scytosiphon-promiscuus_asm6_MALE"
    # "Scytosiphon-promiscuus_mrf5_FEMALE"
    # "Undaria-pinnatifida_phaeo216_FEMALE"
    # "Undaria-pinnatifida_phaeo218_MALE"
    "Ectocarpus-sp7_Ec25_FEMALE"
    # "Ectocarpus-sp7_Ec32_MALE"
)

SPECIES=(
    # "Desmarestia-dudresnayi"
    # "Desmarestia-herbacea"
    # "Desmarestia-herbacea"
    # "Schizocladia-ischiensis"
    # "Scytosiphon-promiscuus"
    # "Scytosiphon-promiscuus"
    # "Undaria-pinnatifida"
    # "Undaria-pinnatifida"
    # "Ectocarpus-sp7"
    "Ectocarpus-sp7"
)
MARKS=(
    # "H3K36me3" 
"H3K4me3" 
# "H3K79me2" "H3K9ac" "H4K20me3"
)

# Ensure Conda environment activation works
source ~/.bashrc  # or ~/.bash_profile if necessary
conda activate deeptools || { echo "Failed to activate deeptools"; exit 1; }

# Function to compute matrix and plot profile
process_sample() {
    local sample="$1"
    local species="$2"
    local mark="$3"

    # Define directories
    BED_DIR="$BASE_DIR/gtf-gff-bed-files/bed_files_for_plotProfile_quantiles/${sample}"
    if [[ ! -d "$BED_DIR" ]]; then
        echo "Error: BED directory $BED_DIR does not exist. Skipping..."
        return
    fi

    # Collect BED files (handle missing files correctly)
    shopt -s nullglob
    BED_FILES=($BED_DIR/*.bed)
    if [[ ${#BED_FILES[@]} -eq 0 ]]; then
        echo "Error: No .bed files found for $sample ($species). Skipping..."
        return
    fi

    # Find bigwig files based on species
    if [[ $species == "Ectocarpus-sp7" ]]; then
        BW_FILES=($PROJECT_DIR/${sample}/05_bigwig_bedgraph_TO_USE/*${mark}*merged.log2ratio.bw)
    else
        BW_FILES=($PROJECT_DIR/${sample}/05_samples_merge-rep_peak-call_bigwigs/*${mark}*merged.log2ratio.bw)
    fi

    if [[ ${#BW_FILES[@]} -eq 0 ]]; then
        echo "Error: No .bw files found for $sample ($species). Skipping..."
        return
    fi

    echo "Processing $sample ($species)..."
    echo "Using BED files: ${BED_FILES[@]}"
    echo "Using BW files: ${BW_FILES[@]}"

    # Compute matrix
    MATRIX_FILE="${OUTPUT_DIR}/matrix_quantiles_${sample}_${mark}.gz"
    computeMatrix scale-regions -S "${BW_FILES[@]}" \
        -R "${BED_FILES[@]}" \
        --beforeRegionStartLength 1000 \
        --afterRegionStartLength 1000 \
        --regionBodyLength 5000 \
        -bs 100 \
        --skipZeros \
        -o "$MATRIX_FILE" || { echo "computeMatrix failed, skipping..."; return; }

    # Plot profile with quantile labels
    PLOT_FILE="${OUTPUT_DIR}/plotProfile_quantiles_${sample}_${mark}.pdf"
    plotProfile -m "$MATRIX_FILE" \
        -out "$PLOT_FILE" \
        --startLabel="start" \
        --endLabel="end" \
        --yAxisLabel="log2(IP/H3)" \
        --plotTitle "Marks on genes in ${sample}" \
        --plotType se \
        --colors "#f8a68d" "#f36f57" "#e43a30" "#a31d21" grey
        # --samplesLabel "1st quantile" "2nd quantile" "3rd quantile" "4th quantile" "no expression"

}

# Loop over samples and process each one
for i in "${!SAMPLES[@]}"; do
    for j in "${!MARKS[@]}"; do
        process_sample "${SAMPLES[$i]}" "${SPECIES[$i]}" "${MARKS[$j]}"
    done
done

# Deactivate conda environment
conda deactivate

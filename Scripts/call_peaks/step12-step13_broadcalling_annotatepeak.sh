#!/bin/bash
# Reserve CPUs for this job
#$ -pe parallel 16
# Request RAM
#$ -l h_vmem=64G
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#$ -o peakcall-merge-annotate.out
#Run job from current working directory
#$ -cwd
#
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas
# Reset the annotatePeaks.out file before each run
> peakcall-merge-annotate-verbose.out

#####################################################################
### Step 12: Recalling with broad option                          ###
### Step 13: Peak annotation                                      ###
#####################################################################

# Call peaks on the merged BAM files
call_peaks() {
  export TMPDIR="/ebio/scratch/jvigneau/deeptool_tmp/"
  local species="$1"
  local species2="$2"
  local output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/${species}/12_broadPeakcalling"
  mkdir -p "$output"
  local input="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/${species}/05_samples_merge-rep_peak-call_bigwigs"
  conda activate macs2
  echo "Step 12: Call peaks on the merged BAM files for ${species} with broadPeak"
  for merged_bam in "$input"/*K*_merged.bam; do
    NAME=$(basename "$merged_bam" _merged.bam)
    # Call all samples with the "broad parameter"
    # macs2 callpeak -n "${NAME}_merged-peaks" -t "$merged_bam" -c "${input}/${species2}_H3_merged.bam" -f BAMPE --outdir "$output" --bdg --broad &
    macs2 callpeak -n "${NAME}_merged-peaks" -t "$merged_bam" -c "${input}/${species2}_H3_merged.bam" -f BAMPE --outdir "$output" --broad &
  done
  wait
  if [ $? -eq 0 ]; then
    echo "Peaks called for ${species}."
  else
    echo "Error occurred during peak calling for ${species}."
  fi
  conda deactivate
}

# Function to annotate peaks
annotate_peaks() {
  local species="$1"
  local species2="$2"
  local genome_assembly="$3"
  local gff_file="$4"
  local peaks_dir="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/${species}/12_broadPeakcalling"
  local output_dir="${projectDir}/${species}/13_annotatePeaks"
  mkdir -p "$output_dir"
  conda activate homer_venv
  for peaks_file in "$peaks_dir"/*_merged-peaks_peaks.broadPeak; do
    local output_file="${output_dir}/$(basename ${peaks_file%_peaks.broadPeak})_annotatePeaks.txt"
    echo "Annotating peaks for file: $peaks_file"
    # Annotate peaks with gene information using HOMER
    annotatePeaks.pl "$peaks_file" "$genome_assembly" -gtf "$gff_file" > "$output_file" &
    # annotatePeaks.pl "$peaks_file" "$genome_assembly" -gid -gff "$gff_file" > "$output_file" &
    # annotatePeaks.pl $peaks_file $genome_assembly -gid -gff $gff_file -gene > $output_file &

    # Check if the annotation was successful
    if [ $? -eq 0 ]; then
      echo "Annotation completed successfully for $peaks_file. Output saved to $output_file"
    else
      echo "Error occurred during annotation for $peaks_file."
    fi
    wait
  done
  conda deactivate
}

### SCRIPT ###
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo"  # Replace with your actual project directory
# Ensure the project directory is available to the script
export projectDir=$PROJECT_DIR



call_peaks "Desmarestia-herbacea_SxS46_FEMALE" "Desmarestia_herbacea_female_SxS46"
annotate_peaks "Desmarestia-herbacea_SxS46_FEMALE" \
  "Desmarestia_herbacea_female_SxS46" \
  "/ebio/abt5_projects/reference_sequences/Desmarestia/Desmarestia_herbacea/genome/reference/Desmarestia-herbacea.reference.fasta" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Desmarestia-herbacea_SxS46_FEMALE/01_results_nf-core-chipseq/genome/Desmarestia-herbacea.reference.gt.gtf"

call_peaks "Desmarestia-dudresnayi_SxS43" "Desmarestia_dudresnayi_SxS43"
annotate_peaks "Desmarestia-dudresnayi_SxS43" \
  "Desmarestia_dudresnayi_SxS43" \
  "/ebio/abt5_projects/reference_sequences/Desmarestia/Desmarestia_dudresnayi/genome/reference/Desmarestia-dudresnayi.chr.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/gtf-gff-bed-files/Desmarestia-dudresnayi.cleaned2025.jv.gffread.gtf" &

call_peaks "Desmarestia-herbacea_SxS47_MALE" "Desmarestia_herbacea_male_SxS47"
annotate_peaks "Desmarestia-herbacea_SxS47_MALE" \
  "Desmarestia_herbacea_male_SxS47" \
  "/ebio/abt5_projects/reference_sequences/Desmarestia/Desmarestia_herbacea/genome/reference/Desmarestia-herbacea.reference.fasta" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Desmarestia-herbacea_SxS46_FEMALE/01_results_nf-core-chipseq/genome/Desmarestia-herbacea.reference.gt.gtf" &

call_peaks "Ectocarpus-sp7_ORO_FEMALE" "OROF"
annotate_peaks "Ectocarpus-sp7_ORO_FEMALE" \
  "OROF" \
  "/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/genome/alternative_versions/Ec32_V5.final.bac_filtered.liniar_T2T.reference.fasta" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_ORO_FEMALE/01_results_nf-core-chipseq/genome/Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gtf" &

call_peaks "Ectocarpus-sp7_ORO_MALE" "OROM"
annotate_peaks "Ectocarpus-sp7_ORO_MALE" \
  "OROM" \
  "/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/genome/alternative_versions/Ec32_V5.final.bac_filtered.liniar_T2T.reference.fasta" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_ORO_FEMALE/01_results_nf-core-chipseq/genome/Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gtf" &

call_peaks "Scytosiphon-promiscuus_mrf5_FEMALE" "Scytosiphon_promiscuus_female_mr5f"
annotate_peaks "Scytosiphon-promiscuus_mrf5_FEMALE" \
  "Scytosiphon_promiscuus_female_mr5f" \
  "/ebio/abt5_projects/reference_sequences/Scytosiphon/Scytosiphon_promiscuus/genome/reference/Scytosiphon-promiscuus_MALE.chr.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/gtf-gff-bed-files/Scytosiphon-promiscuus_MALE_chr_2_AGAT.gtf" &

call_peaks "Scytosiphon-promiscuus_asm6_MALE" "Scytosiphon_promiscuus_male_Asm6"
annotate_peaks "Scytosiphon-promiscuus_asm6_MALE" \
  "Scytosiphon_promiscuus_male_Asm6" \
  "/ebio/abt5_projects/reference_sequences/Scytosiphon/Scytosiphon_promiscuus/genome/reference/Scytosiphon-promiscuus_MALE.chr.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/gtf-gff-bed-files/Scytosiphon-promiscuus_MALE_chr_2_AGAT.gtf" &

call_peaks "Schizocladia-ischiensis_KU333" "Schizocladia-ischiensis-KU333"
annotate_peaks "Schizocladia-ischiensis_KU333" \
  "Schizocladia-ischiensis-KU333" \
  "/ebio/abt5_projects/reference_sequences/Schizocladia/Schizocladia_ischiensis/genome/reference/Schizocladia-ischiensis.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/gtf-gff-bed-files/Schizocladia-ischiensis.gtf" &

call_peaks "Undaria-pinnatifida_phaeo216_FEMALE" "Undaria_pinnatifida_female_PHAEO216"
annotate_peaks "Undaria-pinnatifida_phaeo216_FEMALE" \
  "Undaria_pinnatifida_female_PHAEO216" \
  "/ebio/abt5_projects/reference_sequences/Undaria/Undaria_pinnatifida/genome/reference/PUBLIC_Undaria-pinnatifida_MALE.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/gtf-gff-bed-files/PUBLIC_Undaria-pinnatifida_MALE.jv3.gtf" &

call_peaks "Undaria-pinnatifida_phaeo218_MALE" "Undaria_pinnatifida_male_PHAEO218"
annotate_peaks "Undaria-pinnatifida_phaeo218_MALE" \
  "Undaria_pinnatifida_male_PHAEO218" \
  "/ebio/abt5_projects/reference_sequences/Undaria/Undaria_pinnatifida/genome/reference/PUBLIC_Undaria-pinnatifida_MALE.fa" \
  "/ebio/abt5_projects/PhaeoChromo/data/gtf-gff-bed-files/PUBLIC_Undaria-pinnatifida_MALE.jv3.gtf"
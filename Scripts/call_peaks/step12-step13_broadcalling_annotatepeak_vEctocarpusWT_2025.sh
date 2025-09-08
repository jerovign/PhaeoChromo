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
set -x

#####################################################################
### Step 12: Recalling with broad option                          ###
### Step 13: Peak annotation                                      ###
#####################################################################

genome_assembly="/ebio/abt5_projects/reference_sequences/Ectocarpus/Ec32/genome/reference/Ec32_V5.final.bac_filtered.liniar_T2T.reference.fasta"
gff_file="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_ORO_FEMALE/01_results_nf-core-chipseq/genome/Ec32_V5.final.bac_filtered.liniar_T2T.chr_scaf_sdr.reference.gt.gtf"

#############################################
####################H3K72ME2#################
#############################################

# ### "Ectocarpus-sp7_Ec25_FEMALE"
# PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/05_samples_merge-Rep_peak-call_bigwigs_michael_paired-end/"  # Replace with your actual project directory
# output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/13_annotatePeaks"
# mkdir -p $output
# # Call all samples with the "broad parameter"
# conda activate macs2
#   macs2 callpeak -n "Ectocarpus-sp7_Ec25_FEMALE_H3K79me2_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3K79me2_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3_merged.bam" -f BAMPE --outdir "$output" --broad
# conda deactivate
# # Function to annotate peaks
# conda activate homer_venv
#   annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec25_FEMALE_H3K79me2_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec25_FEMALE_H3K79me2_merged-peaks_annotatePeaks_gene.txt"
# conda deactivate

# ### "Ectocarpus-sp7_Ec32_MALE"
# PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/05_samples_merge-Rep_peak-call_bigwigs_michael_paired-end/"  # Replace with your actual project directory
# output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/13_annotatePeaks"
# mkdir -p $output
# # Call all samples with the "broad parameter"
# conda activate macs2
#   macs2 callpeak -n "Ectocarpus-sp7_Ec32_MALE_H3K79me2_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3K79me2_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3_merged.bam" -f BAMPE --outdir "$output" --broad
# conda deactivate
# # Function to annotate peaks
# conda activate homer_venv
#   annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec32_MALE_H3K79me2_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec32_MALE_H3K79me2_merged-peaks_annotatePeaks_gene.txt"
# conda deactivate

#############################################
####################H3K4ME3#################
#############################################

### "Ectocarpus-sp7_Ec25_FEMALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/05_samples_merge-Rep_peak-call_bigwigs_michael_paired-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec25_FEMALE_H3K4me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3K4me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3_merged.bam" -f BAMPE --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec25_FEMALE_H3K4me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec25_FEMALE_H3K4me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

### "Ectocarpus-sp7_Ec32_MALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/05_samples_merge-Rep_peak-call_bigwigs_michael_paired-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec32_MALE_H3K4me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3K4me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3_merged.bam" -f BAMPE --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec32_MALE_H3K4me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec32_MALE_H3K4me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

#############################################
####################H3K9ac###################
#############################################

### "Ectocarpus-sp7_Ec25_FEMALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec25_FEMALE_H3K9ac_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3K9ac_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec25_FEMALE_H3K9ac_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec25_FEMALE_H3K9ac_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

### "Ectocarpus-sp7_Ec32_MALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec32_MALE_H3K9ac_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3K9ac_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec32_MALE_H3K9ac_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec32_MALE_H3K9ac_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

#############################################
####################H3K36me3#################
#############################################

### "Ectocarpus-sp7_Ec25_FEMALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec25_FEMALE_H3K36me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3K36me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec25_FEMALE_H3K36me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec25_FEMALE_H3K36me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

### "Ectocarpus-sp7_Ec32_MALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec32_MALE_H3K36me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3K36me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec32_MALE_H3K36me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec32_MALE_H3K36me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

#############################################
####################H4K20me3#################
#############################################

### "Ectocarpus-sp7_Ec25_FEMALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec25_FEMALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec25_FEMALE_H4K20me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H4K20me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec25_FEMALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec25_FEMALE_H4K20me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec25_FEMALE_H4K20me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

### "Ectocarpus-sp7_Ec32_MALE"
PROJECT_DIR="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/05_samples_merge-Rep_peak-call_bigwigs_josselin_single-end/"  # Replace with your actual project directory
output="/ebio/abt5_projects/PhaeoChromo/data/early-stage_evochromo/Ectocarpus-sp7_Ec32_MALE/13_annotatePeaks"
mkdir -p $output
# Call all samples with the "broad parameter"
conda activate macs2
  macs2 callpeak -n "Ectocarpus-sp7_Ec32_MALE_H4K20me3_merged-peaks" -t "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H4K20me3_merged.bam" -c "${PROJECT_DIR}/Ectocarpus-sp7_Ec32_MALE_H3_merged.bam" --outdir "$output" --broad
conda deactivate
# Function to annotate peaks
conda activate homer_venv
  annotatePeaks.pl "${output}/Ectocarpus-sp7_Ec32_MALE_H4K20me3_merged-peaks_peaks.broadPeak" "$genome_assembly" -gid -gtf "$gff_file" > "/ebio/abt5_projects/PhaeoChromo/data/box_plots/output/Ectocarpus-sp7_Ec32_MALE_H4K20me3_merged-peaks_annotatePeaks_gene.txt"
conda deactivate

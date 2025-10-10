# Data Preprocessing: Historical *Pseudomonas viridiflava* ATUE5 Genomes

This pipeline describes the steps used to extract, authenticate, and analyze historical *Pseudomonas viridiflava* ATUE5 genomes from *Arabidopsis thaliana* herbarium-derived metagenomes collected between 1817 and 2015.

---

## [Step 1: Extract ATUE5-Mapped Reads](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f07c60d0ca494ab77376b9d2cd0f6163e3f5a068/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h40genomes_stats)

A total of 400+ globally distributed *A. thaliana* herbarium specimens were screened for the presence of *P. viridiflava* ATUE5.

1. **Host depletion**:
   - Raw metagenomic reads were mapped to the *A. thaliana* TAIR10 reference genome using BWA aln.
   - Unmapped reads were retained.

2. **ATUE5 mapping**:
   - Host-depleted reads were mapped to the *P. viridiflava* ATUE5 reference genome (p25.C2).
   - Samples with ≥63% genome breadth at sufficient depth (>1.2x) were retained.
   - **Output**: [46 high-quality historical *Pseudomonas* sp. genomes](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/75ffb2bde41f0ec50d9b0f81b5d4ecd85b16719d/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh46_metadata_withdedup.txt). 

---

## [Step 2: Ancient DNA Authentication](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f07c60d0ca494ab77376b9d2cd0f6163e3f5a068/scripts/step0_datapreprocessing_tailocin_presence/step2_mapdamage)

Ancient DNA authenticity was assessed using mapDamage:

- Quantified 5′ C→T and 3′ G→A substitution frequencies.
- Verified fragment length distributions.
- All 46 samples showed characteristic ancient DNA patterns. [See DNA damage results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step3_h46ps_at), [See first base C to T freq correlation between two speices](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step4_correlation),  [See length distributions](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize). 

---

## [Step 3: Historical ATUE5 phlogeny reconstruction](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f07c60d0ca494ab77376b9d2cd0f6163e3f5a068/scripts/step0_datapreprocessing_tailocin_presence/step3_h46_m85_phylogeny_reconstruction)

To confirm strain identity:

- A maximum likelihood phylogeny was built using 296 biallelic SNPs from the 46 genomes, alongside 55 modern ATUE5 and 30 non-ATUE5 reference genomes.
- Tree inference was performed with IQ-TREE using the TPM3+ASC+R2 substitution model.
- Samples clustering within the ATUE5 clade were classified as ATUE5. [See phylogeny results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step3_h46_m85_tree).
- **Output**: [40 authenticated historical ATUE5 genomes](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/75ffb2bde41f0ec50d9b0f81b5d4ecd85b16719d/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh40_metadata_withdedup.txt).  

---

## [Step 4: Tailocin Region Detection](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f07c60d0ca494ab77376b9d2cd0f6163e3f5a068/scripts/step0_datapreprocessing_tailocin_presence/step4_tailocin_present_stats.sh)

We assessed the presence of the tailocin gene cluster:

- All 40 authenticated genomes showed strong coverage of the tailocin region. [See all details](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/a4169a5d0b2bcdaa69054731d4433f80f4b8e661/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence/tailocin_coverage_summary.tsv).
  - **Average covered proportion**: 0.81  
  - **Average depth**: 28.57×  


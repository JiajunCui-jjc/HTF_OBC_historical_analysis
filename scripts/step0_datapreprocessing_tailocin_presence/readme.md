# Data Preprocessing: Historical *Pseudomonas viridiflava* ATUE5 Genomes

This pipeline describes the steps used to extract, authenticate, and analyze historical *Pseudomonas viridiflava* ATUE5 genomes from *Arabidopsis thaliana* herbarium-derived metagenomes collected between 1817 and 2015.

---

## [Step 1: Extract ATUE5-Mapped Reads](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/629947e97eba44e22798f327b9265f724996d58f/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h49_m57genomes_stats)

A total of 49 globally distributed *A. thaliana* herbarium specimens were screened for the presence of *P. viridiflava* ATUE5.

1. **Host depletion**:
   - Raw metagenomic reads were mapped to the *A. thaliana* TAIR10 reference genome using BWA aln.
   - Unmapped reads were retained.

2. **ATUE5 mapping**:
   - Host-depleted reads were mapped to the *P. viridiflava* ATUE5 reference genome (p25.C2).
   - Samples with ≥60% genome breadth at sufficient depth (>1x) were retained.
   - **Output**: [49 high-quality historical *Pseudomonas* sp. genomes](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/de9d7ced1c979795f5860cd4d0868c7ade7236d0/results/step0_datapreprocessing_tailocin_presence/step1_h46metadata/SuppTable_h46_dedup_updated.txt). 

---

## [Step 2: Ancient DNA Authentication](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/629947e97eba44e22798f327b9265f724996d58f/scripts/step0_datapreprocessing_tailocin_presence/step2_mapdamage)

Ancient DNA authenticity was assessed using mapDamage:

- Quantified 5′ C→T and 3′ G→A substitution frequencies.
- Verified fragment length distributions.
- All 49 samples showed characteristic ancient DNA patterns. [See DNA damage results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step3_h46ps_at), [See first base C to T freq correlation between two speices](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step4_correlation),  [See length distributions](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize). 

---

## [Step 3: Historical ATUE5 phlogeny reconstruction](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/629947e97eba44e22798f327b9265f724996d58f/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction)

To confirm strain identity:

- A maximum likelihood phylogeny was built using 3357 biallelic SNPs (allowed 5% missing info) from the 49 historical genomes, alongside 53 modern ATUE5 and 30 non-ATUE5 reference genomes.
- Tree inference was performed with IQ-TREE using the TVM+F+ASC+R3 substitution model.
- 43 historical isolates clustered within the modern ATUE5 clade were classified as ATUE5. [See phylogeny results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9ba8ed4eff3cf87279332d113c72825d0778153f/results/step0_datapreprocessing_tailocin_presence/step3_h46_m85_tree).
- A maximum likelihood phylogeny was built using 8228 biallelic SNPs (full info sites that appear across all isolates) from the 43 historical ATUE5 genomes and 53 modern ATUE5 genomes.
- **Output**: [43 authenticated historical ATUE5 genomes tree](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/75ffb2bde41f0ec50d9b0f81b5d4ecd85b16719d/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh40_metadata_withdedup.txt).  [43 authenticated historical ATUE5 genomes metadata](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/75ffb2bde41f0ec50d9b0f81b5d4ecd85b16719d/results/step0_datapreprocessing_tailocin_presence/step1_h40metadata/allh40_metadata_withdedup.txt).  

---

## [Step 4: Tailocin Region Detection](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/629947e97eba44e22798f327b9265f724996d58f/scripts/step0_datapreprocessing_tailocin_presence/step4_tailocin_present_stats.sh)

We assessed the presence of the tailocin gene cluster:

- All 43 authenticated genomes showed strong coverage of the tailocin region. [See all details](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/f50945e6ddef9a45a48ec6979ea13be1f51d7572/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence/tailocin_coverage_summary.tsv).
  - **Average covered proportion**: 0.856  
  - **Average depth**: 17.3×  


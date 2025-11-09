# Data Preprocessing: Historical *Pseudomonas viridiflava* ATUE5 Genomes

This pipeline describes the steps used to extract, authenticate, and analyze historical *Pseudomonas viridiflava* ATUE5 genomes from *Arabidopsis thaliana* herbarium-derived metagenomes collected between 1817 and 2015.

---
## [Step 1: Extract ATUE5-Mapped Reads](/scripts/step0_datapreprocessing_tailocin_presence/step1_extract_h49_m57genomes_stats)

A total of 49 globally distributed *A. thaliana* herbarium specimens were screened for the presence of *P. viridiflava* ATUE5.

1. **Host depletion**:
   - Raw metagenomic reads were mapped to the *A. thaliana* TAIR10 reference genome using BWA aln.
   - Unmapped reads were retained.

2. **ATUE5 mapping**:
   - Host-depleted reads were mapped to the *P. viridiflava* ATUE5 reference genome (p25.C2).
   - Samples with ≥60% genome breadth at sufficient depth (>1x) were retained.
   - **Output**: [49 high-quality historical *Pseudomonas* sp. genomes](/results/figureandtable/supptables/SuppTable_h49.txt). 

---

## [Step 2: Ancient DNA Authentication](/scripts/step0_datapreprocessing_tailocin_presence/step2_mapdamage)

Ancient DNA authenticity was assessed using mapDamage:

- Quantified 5′ C→T and 3′ G→A substitution frequencies.
- Verified fragment length distributions.
- All 49 samples showed characteristic ancient DNA patterns. [See DNA damage results](/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot), [See first base C to T freq correlation between two speices](/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/suppfig3mapdamage/49_CtoTcor_Pseudomonas_vs_Arabidopsis.pdf),  [See length distributions](/results/step0_datapreprocessing_tailocin_presence/step2_mapdamageplot/step5_lengthdistribution_insertsize). 

---

## [Step 3: Historical ATUE5 phlogeny reconstruction](/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction)

To confirm strain identity:

- A maximum likelihood phylogeny was built using 3095 biallelic SNPs (allowed 5% missing info) from the 49 historical genomes, alongside 53 modern ATUE5 and 30 non-ATUE5 reference genomes.
- Tree inference was performed with IQ-TREE using the TVM+F+ASC+R4 substitution model.
- 43 historical isolates clustered within the modern ATUE5 clade were classified as ATUE5. [See phylogeny results](/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree).
- A maximum likelihood phylogeny was built using 8228 biallelic SNPs (full info sites that appear across all isolates) from the 43 historical ATUE5 genomes and 53 modern ATUE5 genomes.
- **Output**: [43 authenticated historical ATUE5 genomes tree](/results/step0_datapreprocessing_tailocin_presence/step3_phylogeny_tree/tree_ATUE5).  [43 authenticated historical ATUE5 genomes metadata](/data/historicalATUE5_43.txt).  

---

## [Step 4: Tailocin Region Detection](/scripts/step0_datapreprocessing_tailocin_presence/step4_tailocin_present_stats.sh)

We assessed the presence of the tailocin gene cluster:

- All 43 authenticated genomes showed strong coverage of the tailocin region. [See all details](/results/step0_datapreprocessing_tailocin_presence/step4_tailocin_presence/tailocin_coverage_summary.tsv).
  - **Average covered proportion**: 0.86  
  - **Average depth**: 17.4×  


## [Step 5: LD decay](/scripts/step0_datapreprocessing_tailocin_presence/step3_phylogeny_reconstruction/step1_runh49_m53_30_tree/step4_m53_LDdecay) and [pairwise recombination analysis](/scripts/step0_datapreprocessing_tailocin_presence/step5_recophy_runpairwisevisual.R)

- Genome-wide linkage disequilibrium (LD) and recombination analyses were conducted on 53 modern ATUE5 genomes to test whether HTF–OBC associations arise from physical linkage or selective constraint.
- LD decay was calculated using VCFtools v0.1.17 from 298,977 high-quality biallelic SNPs; pairwise r² values were binned by inter-SNP distance, showing rapid decay within 2 kb.
- Recombination rates were estimated using Recophy, which identified recombinant and clonal SNPs across all pairwise genome comparisons.

- 

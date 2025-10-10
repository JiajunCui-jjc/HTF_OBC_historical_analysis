## Oantigen OBC cluster gene P/A analysis

## historical46_modern57_bybamcov/

From **step1 to step7**, we assessed the presence/absence of six tailocin-associated genes across 40 historical and 57 modern isolates based on BAM coverage.

- In **step1** and **step2**, we identified the coordinates of the six genes within the *Pseudomonas viridiflava* p25.C2 reference genome (which includes multiple contigs).
- In **steps 3–6**, we calculated:
  - **Covered proportion** for each gene
  - **Relative read depth** of each isolate’s mapped library over the gene region
  [See summary table](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_cov_depth_wide_cat.tsv)
> The BAM files were previously generated using `bwa aln`. Example commands are provided in `step0_datapreprocessing`.

- In **step7**, we filtered results using the following thresholds:
  - (i) **Coverage ≥ 50%**
  - (ii) **Mean depth ≥ 75% of the isolate’s genome-wide average**

Genes not meeting both thresholds were scored as **absent**. [The final result is a **binary gene presence/absence matrix** across isolates.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/historical46_modern57_bybamcov/final_sixgenes_binary_matrix.tsv)

---

## checkespE2_contigs/

The **espE2** locus required additional attention due to its size (~4.5 kb) and high sequence divergence.

- For [**modern isolates**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m):
  - Reads mapping to the p25.C2 *espE2* reference were extracted.
  - If a **fully mapped variant to p25.C2 espE2** was found, it was retained.
  - For **partially mapped** hits, the region was extended from the mapped anchor to recover full-length gene sequences (~4,527 ± 200 bp) within the FASTA contigs. Each sequence was manually viewed using Expasy translate to identify start and stop codons. Then sequences were trimmed at the stop codon and then batch-translated into amino acids using EMBOSS Transeq. [modern39_espE2_aa_MSA](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/39fasta_emboss_transeq-I20250724-161752-0560-416209-p1m.out.txt)
  - This yielded a diverse set of modern *espE2* haplotypes.[See msa with 18 empty and 39 retrieved seqs.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt)

- For [**historical isolates**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/fish_historical) (which lack assemblies):
  - The above [**modern haplotype set**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt) was used as an expanded reference.
  - Reads mapping to any haplotype were extracted and **locally assembled** with **SPAdes**.
  - Contigs were then aligned back to the haplotype reference with **Minimap2** to identify the best-matching allele. [See summary binary table including espE2_rescued here](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h36_epsE2rescued_binary_matrix.tsv)

- For **both modern and historical** datasets:
  - *espE2* sequences were **translated into amino acids** using Expasy.
  - **ORFs were identified** via start/stop codon detection.
  - Sequences were **aligned** with **Clustal Omega**.
  - Final multiple sequence alignments and allele length distributions are provided [here](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b73cbe487b45446d374f85cb44b5df5bc910dfac/results/step3_combine/combined_HTF_oantigen_m57_h35.txt).

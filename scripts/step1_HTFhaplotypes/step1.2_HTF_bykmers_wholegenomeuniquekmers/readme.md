### HTF k-mer identification and additional analysis pipeline

---

### [`step0_wholegenome_7refs/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/99001840be68251201d8ad9982588e861863dc92/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step0_wholegenome_7refs)

- `step1_generatepaf.sh`  
- `step2_summarypaf_corrected_byhand.sh`  
- `step3_searchanddelete.sh`  

> In **step0**, I run `minimap2` to find the coordinates for the HTF region in each reference genome, then delete them to generate 7 whole-genome backgrounds that exclude the self HTF region for downstream k-mer calling.

---

### [`step1.1_refkmers_all_to_wgunique.sh`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.1_refkmers_all_to_wgunique.sh)

> In **step1.1**, I run the k-mer calling to generate 7 sets of kmers that are unique from each other and also exclusive from the 7 representative OTU5 whole genome backgrounds.

---

### [`step1.2_refkmers_hamming_filtering/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_refkmers_hamming_filtering)

> In **step1.2**, I compute the hamming distance matrix and visualize it, then use R to iteratively filter out any k-mer pairs with distance ≤1.  
> This yields the final `kmers_unique_hamming2` set. [See results.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/99001840be68251201d8ad9982588e861863dc92/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2) [See Rplot of kmers locations.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/99001840be68251201d8ad9982588e861863dc92/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step1_cal_hamming_pairwise_matrix/Rvisualization/runlocationofkmers/HTF_kmer_positions_after.pdf)

---

### [`step2.1_generateisolate_jf_kmersm57.sh` ](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step2.1_generateisolate_jf_kmersm57.sh)
### [`step2.2_bamtofastq_kmersh49.sh`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step2.2_bamtofastq_kmersh49.sh)
### [`step2.3_generateisolate_jf_kmersh49_arraysub.sh`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step2.3_generateisolate_jf_kmersh49_arraysub.sh)


> Here I generate `.jf` files for isolate k-mers.  
> For H49 samples, I remove duplicates and use the `rmdup` FASTQ to generate `.jf` files.

---

### [`step3_wg_background_distribution/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step3_wg_background_distribution)

- `readme`  
- `step1_dump_wgdistribution_Rvisual.R`  
- `step2_calculate_mean_wg_kmercount_asthreshold.sh`  
- `step3_dump_wgdistribution_Rvisual.R`  

> In **step3**, I analyze whole-genome (WG) k-mer distributions.  
> I visualize the distribution in R ( 2 < read <=mean +3*sd to show higher resolution). [See the whole genome kmer distribution](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/8aec893343d313dd5b48790f80e923c219184e99/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/distribution/filtered)

---

### [`step4_query_and_filterby_wgthreshold/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/8aec893343d313dd5b48790f80e923c219184e99/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold)

> In **step4**, I query each isolate’s k-mers and filter them using the average whole-genome kmer depth threshold (mean – 0.5 SD for historical, mean wg kmer depth for modern) after excluding low-depth wg k-mers (≤ 1). 
> Then I:
> - Summarize the raw k-mer counts of each HTF–isolate pair
> - Use R to calculate the proportion of HTF-specific k-mers to the total matched HTF k-mers for that isolate
> - Use this index to rank the best HTF matches  
> - Generate mixture histograms of HTF and WG k-mer distributions
> To do this, I run `step1.sh` to obtain `wg_p25.c2_depthsummary` (long-format k-mer counts), then run `mix.R` to plot the distribution. [See mixed kmer distribution plots](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/c8321aa3b4a9b4bd0ee00a4b03a66d36255e3c59/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix)

---

### [`step5_additional_analysis/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/b7142d1efcbe7022caad783e09d3235af973e87a/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq/authenticate_the2breakingisolates)

> In **step5_additional_analysis**, I:
> - Run coinfection detection
> - Isolates that break the HTF-OBC pattern:
>   - Remove `64.GBR` as a uncertain calls  
>   - Confirm `PL0240` as a confident break
> - [See analysis summary slides](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c8321aa3b4a9b4bd0ee00a4b03a66d36255e3c59/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step2_coinfection/breaking_coninfection/additional_analysis_breakingisolates_coinfection.pptx)
> - Calculate HTF length group frequency differences, plot boxplots and raw data versus time and geography.

---

### Final Combine

#### `combine/step1/`

> In `combine`, I merge the HTF results with O-antigen data [See table](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c8321aa3b4a9b4bd0ee00a4b03a66d36255e3c59/results/step3_combine/combined_HTF_oantigen_m57_h43.txt).  
> This is the [**final main figure**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/3c0c43dabd8abe6f4712b38e43f595d691cb1f51/results/figureandtable/figs/pdfs/mainfig1.png).

#### `../confidence_matrix/`

> In `confidence_matrix`, I combine **local assembly** and **k-mer** results to produce a **confidence matrix of HTF calls**. [See plots](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/3c0c43dabd8abe6f4712b38e43f595d691cb1f51/results/figureandtable/figs/pdfs/othersupps/supp5_confidence)

---

**Done**

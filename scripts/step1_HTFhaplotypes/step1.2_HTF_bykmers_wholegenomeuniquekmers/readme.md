### HTF k-mer identification and additional analysis pipeline

---

### [`step0_wholegenome_7refs/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/015582c21dc87aebe781ac039d76fdd8c83a46d7/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step0_wholegenome_7refs)

- `step1_generatepaf.sh`  
- `step2_summarypaf_corrected_byhand.sh`  
- `step3_searchanddelete.sh`  

> In **step0**, I run `minimap2` to find the coordinates for the HTF region in each reference genome, then delete them to generate 7 whole-genome backgrounds that exclude the self HTF region for downstream k-mer calling.

---

### [`step1.1_refkmers_all_to_wgunique.sh`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/015582c21dc87aebe781ac039d76fdd8c83a46d7/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.1_refkmers_all_to_wgunique.sh)

> In **step1.1**, I run the k-mer calling to generate 7 sets of kmers that are unique from each other and also exclusive from the 7 representative OTU5 whole genome backgrounds.

---

### [`step1.2_refkmers_hamming_filtering/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/015582c21dc87aebe781ac039d76fdd8c83a46d7/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step1.2_refkmers_hamming_filtering)

> In **step1.2**, I compute the hamming distance matrix and visualize it, then use R to iteratively filter out any k-mer pairs with distance ≤1.  
> This yields the final `kmers_unique_hamming2` set. [See results.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/6e2e399e4e7be69c1d048f882571e1418aed81b2/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2) [See Rplot of kmers locations.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/6e2e399e4e7be69c1d048f882571e1418aed81b2/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/kmers_unique_hamming2/runlocationofkmers/HTF_kmer_positions_after.pdf)

---

### [`step2.1_generateisolate_jf_kmersm57.sh` ](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/f1b328393e4d77bc2a5dd24a1d6d7a1e82b1d0ec/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step2.1_generateisolate_jf_kmersm57.sh)
### [`step2.2_generateisolate_jf_kmersh40.sh`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/f1b328393e4d77bc2a5dd24a1d6d7a1e82b1d0ec/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step2.2_generateisolate_jf_kmersh40.sh)

> Here I generate `.jf` files for isolate k-mers.  
> For H40 samples, I remove duplicates and use the `rmdup` FASTQ to generate `.jf` files.

---

### [`step3_wg_background_distribution/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f1b328393e4d77bc2a5dd24a1d6d7a1e82b1d0ec/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step3_wg_background_distribution)

- `readme`  
- `step1_dump_wgdistribution_Rvisual.R`  
- `step2_calculate_mean_wg_kmercount_asthreshold.sh`  
- `step3_dump_wgdistribution_Rvisual.R`  

> In **step3**, I analyze whole-genome (WG) k-mer distributions.  
> I visualize the distribution in R, then compute the mean WG k-mer count (after quality control, 2 < read <=mean +3*sd ) as a threshold for filtering. [See the whole genome kmer distribution](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/9dcc3b3ae058536717f726ee66d104b129c770fa/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step2_isolate_kmers/isolate_jf/distribution/filtered)

---

### [`step4_query_and_filterby_wgthreshold/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/f1b328393e4d77bc2a5dd24a1d6d7a1e82b1d0ec/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step4_query_and_filterby_wgthreshold)

> In **step4**, I query each isolate’s k-mers and filter them using the average WG k-mer threshold.  
> Then I:
> - Summarize the raw k-mer counts of each HTF–isolate pair
> - Use R to calculate the proportion of HTF-specific k-mers to the total matched HTF k-mers for that isolate
> - Use this index to rank the best HTF matches  
> - Generate mixture histograms of HTF and WG k-mer distributions
> To do this, I run `step1.sh` to obtain `wg_p25.c2_depthsummary` (long-format k-mer counts), then run `mix.R` to plot the distribution. [See mixed kmer distribution plots](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/0be4345d26b12a9e1e9d88831ca1f1476a78677e/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/kmer_distribution_mix/distributions)

---

### [`step5_additional_analysis/`](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/8172ca79b8794439f3dfb40b923af94470a24986/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/step5_additional_analysis_coinfection_Oantigenbreakingisolates_HTFlengthfreq)

> In **step5_additional_analysis**, I:
> - Run coinfection detection
> - Isolates that break the HTF-OBC pattern:
>   - Remove `64.GBR` as a uncertain calls  
>   - Confirm `PL0240` as a confident break
> - [See analysis summary slides](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/1cf701cfecc3025ad6d43a7c47ae4783d15816d2/results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step4_additionalanalysis/step2_coinfection/summary_breaking_coninfection)
> - Calculate HTF length group frequency differences, plot boxplots and raw data versus time and geography.

---

### Final Combine

#### `combine/step1/`

> In `combine`, I merge the HTF results with O-antigen data [See table](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/8f6838eebf2e91b341e8b2133a905cdca54a4d69/results/figures_tables/supptables/supptable2sum/merged_final_HTF_tailocin.txt).  
> This is the [**final main figure**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c2712d1d74848bd6f7a825a79ea9944951e8d4ed/results/figures_tables/suppfigs/pdfs/mainfig1new.png).

#### `../confidence_matrix/`

> In `confidence_matrix`, I combine **local assembly** and **k-mer** results to produce a **confidence matrix of HTF calls**. [See plots](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/c2712d1d74848bd6f7a825a79ea9944951e8d4ed/results/figures_tables/suppfig1_confidence_matrix_HTF/step2_Rvisual)

---

**Done**

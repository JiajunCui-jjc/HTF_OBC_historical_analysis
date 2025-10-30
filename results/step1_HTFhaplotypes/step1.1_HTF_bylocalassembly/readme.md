## HTF Haplotype Recovery by Local Assembly

### Historical Samples

For historical samples, we extracted reads mapped to the reference tailocin region containing seven HTF haplotypes. Assemblies were performed using these reads, and contigs were aligned to the reference using Minimap2.

- [**Step 1:**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/6283a4b76c4be91fd127ae36e1751565d91d751e/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical40/step1_extractreadsfrom46h) Successfully assembled tailocin regions for **[29 out of 40 historical samples](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings/samples29withcontigs.txt)** using locally mapped reads.

- [**Step 2:**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/dbe9a96fada18c8a308773c4212635b4d7dcb03e/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical40/step2_assembly) Identified HTF haplotypes in **24 historical assemblies**.

- **Optional step 3:** For MSA, after applying a **65% coverage threshold** to the best-matching reference haplotype, **10 historical HTFs** were retained.

- **Optional step 4:** One of the nucleotide sequences (`PL0046`) was found to be incomplete after translation to amino acids and was excluded.  

- **[Final count](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/dbe9a96fada18c8a308773c4212635b4d7dcb03e/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical40/step3_summarycovprop_topselect.sh):  24 historical HTFs retrieved (9/24 high-confidence).** [See results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings/h40_HTFhaplotype_bylocalassembly.txt)

Later, HTF identity was further confirmed via k-mer presence analysis, which recovered **35 historical HTFs** (including some not retained here due to quality).

---

### Modern Samples (n = 57)

For modern isolates, we aligned local assemblies to the HTF reference using Minimap2. [See the summary table.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/modern57paf/modern57_HTFhaplotype_bylocalassembly.txt)

- For each sample, the **top-matching HTF and TFA reference** was selected based on **coverage proportion**.
- Most samples showed **>98% coverage** to their top HTF reference, along with presence of the matching TFA region from the same haplotype.
- Occasionally, the best-matching TFA may not correspond to the same HTF as the top hit.
- Final haplotype calls were validated using the **HTF k-mer confidence matrix** and the prior haplotype assignment in (Talia. B, 2024).

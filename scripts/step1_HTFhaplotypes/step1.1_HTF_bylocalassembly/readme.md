## HTF Haplotype Recovery by Local Assembly

### Historical Samples

For historical samples, we extracted reads mapped to the reference tailocin region containing seven HTF haplotypes. Assemblies were performed using these reads, and contigs were aligned to the reference using Minimap2.

- [**Step 1:**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/a89d1ea3d6714c02d1dba31b0509bd75adfde13f/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step1_extractreadsfrom49h) Successfully assembled tailocin regions for **[29 out of 40 historical samples](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings/samples29withcontigs.txt)** using locally mapped reads.

- [**Step 2:**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/a89d1ea3d6714c02d1dba31b0509bd75adfde13f/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step2_assembly) Identified HTF haplotypes in **24 historical assemblies**.

- **note** One of the nucleotide sequences (`PL0046`) was found to be incomplete after translation to amino acids and was excluded.  

- **[Final count](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/a89d1ea3d6714c02d1dba31b0509bd75adfde13f/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step3_summarycovprop_topselect.sh):  24 historical HTFs retrieved (9/24 high-confidence).** [See results](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical46paf/mappings/h40_HTFhaplotype_bylocalassembly.txt)

Later, HTF identity was further confirmed via k-mer presence analysis, which recovered **35 historical HTFs** (including some not retained here due to quality).

---

### Modern Samples (n = 57)

[For modern isolates, we aligned local assemblies to the HTF reference using Minimap2](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/497ac25d1406030dd896b68ba2cfc41c1e8f0484/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/modern57/runminimapandselecttopmatched_formodern57.sh). [See the summary table.](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/1d83183b8e188945862609aa0d7e956b8eb24402/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/modern57paf/modern57_HTFhaplotype_bylocalassembly.txt)



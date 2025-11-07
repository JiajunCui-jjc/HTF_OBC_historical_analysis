## HTF Haplotype Recovery by Local Assembly

### Historical Samples

For historical samples, we extracted reads mapped to the reference tailocin region containing seven HTF haplotypes. Assemblies were performed using these reads, and contigs were aligned to the reference using Minimap2.

- [**Step 1:**](/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step1_extractreadsfrom49h) Successfully assembled tailocin regions for **[27 out of 43 historical samples](/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/mappings/sampleswithassemblies.txt)** using locally mapped reads.

- [**Step 2:**](/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step2_assembly) Identified HTF haplotypes in **26 historical assemblies**.

- **[Final count](/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49/step3_summarycovprop_topselect.sh):  26 historical HTFs retrieved.** [See results](/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/mappings/h49_HTFhaplotype_bylocalassembly.txt)

Later, HTF identity was further confirmed via k-mer presence analysis, which recovered **38 historical HTFs**.

---

### Modern Samples (n = 57)

[For modern isolates, we aligned local assemblies to the HTF reference using Minimap2](/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/modern57/runminimapandselecttopmatched_formodern57.sh). [See the summary table.](/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/modern57paf/modern57_HTFhaplotype_bylocalassembly.txt)



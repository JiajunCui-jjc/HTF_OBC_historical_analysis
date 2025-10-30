## HTF–O–Antigen–Gene Presence Summary and espE2 recovering

## [Step 1:](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/3d6f964db5d50e89721a8571a4ffec01c5bb03a5/scripts/step3_combine/step1_generateasummarytable_HTF_Oantigen.sh)

## Step 1.1: Input Files

### (a) OBC-Gene Binary Matrix (including espE2 and espE2 length)

**File path:**

[./results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h49_epsE2rescued_binary_matrix.tsv](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c46237ff992b723c996337e17c086a52cd79cd09/results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h49_epsE2rescued_binary_matrix.tsv)

notes for [espE2 length calculation](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/3d6f964db5d50e89721a8571a4ffec01c5bb03a5/scripts/step3_combine/step1_generateasummarytable_HTF_Oantigen.sh):
- If a sequence is present → espE2_PA = 1  
- If absent → espE2_PA = 0  
- Length is calculated from the FASTA entry  

### (b) HTF Haplotypes Table

**File path:**

[./results/step3_combine/HTF_haplotypes/form57_h49_HTF_bykmer.txt](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c46237ff992b723c996337e17c086a52cd79cd09/results/step3_combine/HTF_haplotypes/form57_h49_HTF_bykmer.txt)

Merged by /results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/h49sample_bestHTF.txt and /results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/m57sample_bestHTF.txt

---

## Step 1.2: Derive O-Antigen Classification from HTF

HTF Haplotypes → OBC Classification  
*1830, *1383 → OBC_absent  
*1803, *1245 → OBC_present

Column name: HTFgroup_Oantigen_PA

---

## Step 1.3: Output Summary Table

Combined columns:  
sample, HTFgroup_Oantigen_PA, wfgD_PA, rmlC1_PA, tagG1_PA, tagH1_PA, spsA_PA, espE2_PA, espE2length, note

**Output file:**

[/results/step3_combine/combined_HTF_oantigen_m57_h43.txt](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c46237ff992b723c996337e17c086a52cd79cd09/results/step3_combine/combined_HTF_oantigen_m57_h43.txt)

---

## Step 2: R plot 
[final figure, HTF - OBC](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/c46237ff992b723c996337e17c086a52cd79cd09/results/step3_combine/mainfig_HTFoantigen_m53_h43_38) 

---

## Step 3: [Confidence matrix](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/c46237ff992b723c996337e17c086a52cd79cd09/results/step3_combine/suppfig1_confidence_matrix_HTF)
columns: modern57 and h43
rows: HTF by kmer, HTF by local assembly and HTF modern from previous study



## HTF–O–Antigen–Gene Presence Summary and espE2 recovering

## [Step 1:](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/40b14f7905c53ea71f77e7c36d4c439fb4a744b6/scripts/step3_combine/step1_generateasummarytable_HTF_Oantigen.sh)

## Step 1.1: Input Files

### (a) OBC-Gene Binary Matrix (including espE2 and espE2 length)

**File path:**

[./results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h36_epsE2rescued_binary_matrix.tsv](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/f2e0253fdd33c23b7d3d43982dbe2bd895d1e699/results/step2_Oantigengenes/espE2_rescue/fish_historical/final_sixgenes_withm57_h36_epsE2rescued_binary_matrix.tsv)

notes for [espE2 length calculation](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/40b14f7905c53ea71f77e7c36d4c439fb4a744b6/scripts/step3_combine/step1_generateasummarytable_HTF_Oantigen.sh):
[./results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/f2e0253fdd33c23b7d3d43982dbe2bd895d1e699/results/step2_Oantigengenes/espE2_rescue/checkespE2_57m/final_afterexpasy/completeopenframe.18empty_39fasta.txt)

- If a sequence is present → espE2_PA = 1  
- If absent → espE2_PA = 0  
- Length is calculated from the FASTA entry  
- Column note is added manually

### (b) HTF Haplotypes Table

**File path:**

[./results/step3_combine/HTF_haplotypes/form57_h35_HTF_bykmer.txt](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/40b14f7905c53ea71f77e7c36d4c439fb4a744b6/results/step3_combine/HTF_haplotypes/form57_h35_HTF_bykmer.txt)
Merged by /results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/h40sample_bestHTF.txt and /results/step1_HTFhaplotypes/step1.2_HTF_bykmers/step3_query_kmers/m57sample_bestHTF.txt

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

[/results/step3_combine/combined_HTF_oantigen_m57_h35.txt](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/40b14f7905c53ea71f77e7c36d4c439fb4a744b6/results/step3_combine/combined_HTF_oantigen_m57_h35.txt)

---

## Step 2: R plot 
[final figure, HTF - OBC](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/c2712d1d74848bd6f7a825a79ea9944951e8d4ed/results/figures_tables/suppfigs/pdfs/mainfig1new.png) 

---

## Step 3: [Confidence matrix](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/tree/c2712d1d74848bd6f7a825a79ea9944951e8d4ed/results/figures_tables/suppfig1_confidence_matrix_HTF/step2_Rvisual)
columns: modern57 and h35
rows: HTF by kmer, HTF by local assembly and HTF modern from previous study



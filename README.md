
HTF and O-antigen Profiling of ***Pseudomonas viridiflava*** ATUE5 in historical ***Arabidopsis thaliana*** Metagenomes
==================================================================================

This repository contains scripts, data, and results for profiling the hypothetical tail fiber (HTF) haplotypes and O-antigen biosynthesis gene content of ***P. viridiflava*** isolates from historical herbarium and modern ***Arabidopsis thaliana*** samples.

Main Figure
--------------

![HTF - OBC pattern](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/bebf869a0fac5939cead9665bfe84027418704c9/results/figures_tables/suppfigs/pdfs/mainfig1new.png)
[samples included](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/ed44f258896c7e46bad668f5938bbb8344963698/results/step3_combine/mainfig_HTFoantigen_m53_h34/m53_h34_names.txt)

Analysis pipeline
--------------
[**Data Preprocessing**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/3178cc9bf9c24f7eb25c42d21a8df33f2d1ff3da/scripts/step0_datapreprocessing_tailocin_presence/readme.md)
 
- Extract ATUE5-mapped reads from historical plant herbarium metagenomes
- Authenticate historical DNA damage patterns
- Historical ATUE5 phlogeny reconstruction
- Identify the presence of tailocin region in historical ATUE5 genomes

**HTF Haplotype Assignment**:

- [Local Assembly Approach](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/8a2d34315a1ab4712382a89c205be19589da1d3f/scripts/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/readme.md):
    - Extract reads mapping to HTF/TFA regions
    - Assemble with SPAdes (skip for modern samples, since we have modern assemblies)
    - Assign best haplotype based on covered proportion (minimap2)
- [K-mer Based Approach](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/3de19f83bff4fdf09f96a5715ff62b79b53e3f79/scripts/step1_HTFhaplotypes/step1.2_HTF_bykmers_wholegenomeuniquekmers/readme.md):
    - Build whole-genome-exclusive HTF-unique kmers
    - Apply iterative Hamming ≥ 2 filtering across haplotypes
    - Query isolate .jf files (Jellyfish)
    - Assign dominant haplotype by max (HTF matched kmers/total matched kmers) proportion
    - Detect coinfections if multiple haplotypes exceed threshold
    - HTF length group frequency distribution

[**O-antigen Gene P/A Detection**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/e936a5a7c9dd58f8d5f3515d23d2a76afec5b014/scripts/step2_Oantigengenes/readme.md):

- A gene is considered present if:
    (i) coverage ≥ 50%
    (ii) mean depth ≥ 75% of genome-wide average
- espE4 handled separately via extended mapping and contig rescue

[**Combined Analysis**](https://github.com/JiajunCui-jjc/HTF_OBC_historical_analysis/blob/a695d98148d520428c049db639f5bfdcfd774c5d/scripts/step3_combine/readme.md):

- Merge HTF and OBC profiles
- Output metadata tables and combined heatmaps

Large data including raw fastq, fasta, reference and so on are stored on Zenodo. Link: 

[Directory Structure](https://github.com/CJJ8848/HTF_OBC_historical_analysis/blob/01535b6afe4fab72f44866cbcf305e59b1537f87/structure.txt)
-------------------
data/

    - modern57.txt, 40h_OTU5.txt: lists of modern/historical samples.
    - readme.txt: description of input formats.

scripts/

    - step0_datapreprocessing_tailocin_presence/: extract ATUE5-mapped reads from historical plant herbarium metagenomes, historical DNA authentication, ATUE5 phylogeny reconstruction and identify tailocin regions in ATUE5 genomes
    - step1_HTFhaplotypes/: HTF haplotype detection
        - step1.1_HTF_bylocalassembly/: local assembly and mapping-based assignment
        - step1.2_HTF_bykmers_wholegenomeuniquekmers/: kmer filtering, querying, and assignment
    - step2_Oantigengenes/: detection of six O-antigen genes P/A including espE2 rescue
    - step3_combine/: integration of HTF and O-antigen data, generation of summary and plots

results/

    - figures_tables/: Contains all figures and tables used in the manuscript and all analyses
    - step0_datapreprocessing_tailocin_presence/: h40 metadata, historical DNA authentication and summary table of tailocin presence in ATUE5 genomes
    - step1_HTFhaplotypes/: HTF haplotype results (assembly/kmer-based)
    - step2_Oantigengenes/: binary P/A gene matrix and espE2 analysis
    - step3_combine/: combined tables and plots (for manuscript figures)


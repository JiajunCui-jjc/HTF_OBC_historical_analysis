# Tailocin Genomic Similarity Analysis — README

This README describes the computational workflow used to assess tailocin genomic similarity across *Pseudomonas viridiflava* ATUE5 (n=53) and *Pseudomonas syringae* (n=126) genomes. Visualisations (PCoA plots and NJ trees) are coloured by **tail fiber type** (T1/T2) and **phylogroup**.

---

## 1. *Pseudomonas syringae* Genome Dataset

Publicly available *P. syringae* genomes were downloaded from NCBI (source: Patricia, TBC) and filtered for assembly quality:

- Assembly length: **5.5–7.5 Mb**
- Number of contigs: **< 1,000**
- N50: **> 100 kb**

To reduce redundancy, genomes were clustered at **99% ANI** and one representative per cluster was selected based on best assembly metrics, yielding a final dataset of **126 *P. syringae* genomes**. For each genome, the tailocin gene cluster flanked by the *trpE* and *trpG* islands was extracted based on GenBank annotations and validated with Pharokka (Bouras et al., 2023; https://doi.org/10.1093/bioinformatics/btac776) and PHOLD (Bouras et al., 2026; https://doi.org/10.1093/nar/gkaf1448).

---

## 2. Tailocin Region Extraction

Tailocin region sequences were extracted from 53 modern assemblies. The tailocin region was defined from the end of *trpE* to the start of *trpG*, using the *P25.C2* assembly as reference:

- **Reference query (minimap2):** `utg000001l_p25.C2:2695984-2717365`
- **Query length:** 21,382 bp

Each modern assembly was searched for the best minimap2 match to this reference query. In cases where tailocin genes spanned multiple contigs, the matched segments were concatenated to reconstruct the full tailocin sequence. The resulting sequences were exported as a multi-FASTA for downstream inspection and phylogenetic analyses. HTF sequences were extracted separately and open reading frames were manually inspected for both nucleotide and amino acid sequences.

---

## 3. k-mer Indexing with panKmer

Three separate indices were built from nucleotide sequences of the tailocin gene cluster:

| Dataset | Description |
|---|---|
| `ALL_genes` | Full tailocin gene cluster (including HTF) |
| `without_HTF` | Tailocin gene cluster excluding the tail fiber gene (HTF) |
| `only_HTF` | Tail fiber gene (HTF) alone |

```bash
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/ALL_genes    -o ALL_genes.tar
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/without_HTF  -o without_HTF.tar
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/only_HTF     -o only_HTF.tar
```

> k-mer length: **31 bp**

---

## 4. Adjacency Matrices

Adjacency matrices were computed from each index:

```bash
pankmer adj-matrix -i ALL_genes.tar    -o ALL_genes.adjmatrix.tsv
pankmer adj-matrix -i without_HTF.tar  -o without_HTF.adjmatrix.tsv
pankmer adj-matrix -i only_HTF.tar     -o only_HTF.adjmatrix.tsv
```

A heatmap with Jaccard distances was also generated for visual inspection:

```bash
pankmer clustermap -i ALL_genes.adjmatrix.tsv -o ALL_genes.adjmatrix.jaccard.svg \
    --metric jaccard --width 20 --height 20
```

---

## 5. PCoA (Principal Coordinates Analysis)

Jaccard distances were computed from each adjacency matrix. PCoA was then performed using the `pcoa` function from the **scikit-bio** Python package.

Points in PCoA figures are coloured by:
- **Tail fiber type** (T1 / T2) — see Section 8
- **Phylogroup** — see Section 8

The full analysis code is in the **Jupyter Notebook** (see `notebooks/`).

---

## 6. Neighbour-Joining (NJ) Trees

Unrooted NJ trees were reconstructed from Jaccard distance matrices using the `DistanceTreeConstructor` function from **Biopython**. Trees are visualised in two layouts: **radial** and **polar** (using a Python script and R).

Output tree files:

```
NJ_trees/ALL_genes.distance.jaccard.NJ.newick
NJ_trees/only_HTF.distance.jaccard.NJ.newick
NJ_trees/without_HTF.distance.jaccard.NJ.newick
```

Code for tree construction and visualisation is in the Jupyter Notebook and the accompanying R/Python scripts.

---

## 7. Tail Fiber Type Assignment (T1 / T2)

Tail fibers were classified as **T1** or **T2** based on BLASTp similarity to 10 reference proteins from Fautt et al. (2025).

- Each *P. syringae* tail fiber (query) was searched against all 10 reference proteins.
- Top hits were ranked by **bitscore**.
- Sequences that did not pass an e-value threshold of **10⁻⁵** were excluded.

---

## 8. Phylogroup Assignment

Phylogroups were assigned to *P. syringae* strains using **ANI** comparisons against reference strains from Marques et al. (2024) (https://www.nature.com/articles/s41597-024-03003-x).

- A genome is assigned to a phylogroup if ANI ≥ **95%** with the corresponding reference.
- If ANI ≥ 95% with multiple references, the **highest-scoring** reference is chosen.
- If ANI < 95% with all references, the genome is **unassigned** (best-hit reference retained). ~10 genomes fall into this category.

---

## Software & References

| Tool | Reference |
|---|---|
| minimap2 | Li (2018) |
| Pharokka | Bouras et al. (2023) — https://doi.org/10.1093/bioinformatics/btac776 |
| PHOLD | Bouras et al. (2026) — https://doi.org/10.1093/nar/gkaf1448 |
| panKmer | Aylward et al. (2023) |
| scikit-bio (`pcoa`) | Aton et al. (2026) |
| Biopython (`DistanceTreeConstructor`) | Cock et al. (2009) |
| BLASTp | NCBI BLAST |
| Tail fiber reference proteins | Fautt et al. (2025) |
| Phylogroup reference strains | Marques et al. (2024) |

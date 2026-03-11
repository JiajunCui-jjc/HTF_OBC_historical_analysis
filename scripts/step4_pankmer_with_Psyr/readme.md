# Tailocin Genomic Similarity Analysis — README

This README describes the computational workflow used to assess tailocin genomic similarity across *Pseudomonas viridiflava* ATUE5 (n=53) and *Pseudomonas syringae* (n=126) genomes. Visualisations (PCoA plots and NJ trees) are coloured by **tail fiber type** (T1/T2) and **phylogroup**.

---

## 1. k-mer Indexing with panKmer

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

## 2. Adjacency Matrices

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

## 3. PCoA (Principal Coordinates Analysis)

Jaccard distances were computed from each adjacency matrix. PCoA was then performed using the `pcoa` function from the **scikit-bio** Python package.

Points in PCoA figures are coloured by:
- **Tail fiber type** (T1 / T2) — see Section 5
- **Phylogroup** — see Section 6

The full analysis code is in the **Jupyter Notebook** (see `notebooks/`).

---

## 4. Neighbour-Joining (NJ) Trees

Unrooted NJ trees were reconstructed from Jaccard distance matrices using the `DistanceTreeConstructor` function from **Biopython**. Trees are visualised in two layouts: **radial** and **polar** (using a Python script and R).

Output tree files:

```
NJ_trees/ALL_genes.distance.jaccard.NJ.newick
NJ_trees/only_HTF.distance.jaccard.NJ.newick
NJ_trees/without_HTF.distance.jaccard.NJ.newick
```

Code for tree construction and visualisation is in the Jupyter Notebook and the accompanying R/Python scripts.

---

## 5. Tail Fiber Type Assignment (T1 / T2)

Tail fibers were classified as **T1** or **T2** based on BLASTp similarity to 10 reference proteins from Fautt et al. (2025), which were used in that study to identify tail fibers across the *P. syringae* species complex.

**Protocol:**
- Each *P. syringae* tail fiber (query) was searched against all 10 reference proteins.
- Top hits were ranked by **bitscore**.
- A tail fiber was excluded from classification if it failed passing an e-value threshold of **10⁻⁵**.

---

## 6. Phylogroup Assignment

Phylogroups were assigned to *P. syringae* strains using **Average Nucleotide Identity (ANI)** comparisons against phylogroup reference strains, following the type strain definitions in **Marques et al. (2024)** (https://www.nature.com/articles/s41597-024-03003-x).

**Assignment rules:**
- A genome is assigned to a phylogroup if ANI ≥ **95%** with the corresponding reference strain.
- If ANI ≥ 95% with multiple references, the **highest-scoring** reference is chosen.
- If ANI < 95% with all references, the genome is **unassigned** (best-hit reference is retained for information). ~10 genomes fall into this category.

---

## Software & References

| Tool | Reference |
|---|---|
| panKmer | Aylward et al. (2023) |
| scikit-bio (`pcoa`) | Aton et al. (2026) |
| Biopython (`DistanceTreeConstructor`) | Cock et al. (2009) |
| BLASTp | NCBI BLAST |
| Tail fiber reference proteins | Fautt et al. (2025) |
| Phylogroup reference strains | Marques et al. (2024) |

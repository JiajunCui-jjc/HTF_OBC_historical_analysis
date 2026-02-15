# 1. panKmer
## 1.1. Create indices
```bash
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/ALL_genes -o ALL_genes.tar
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/without_HTF -o without_HTF.tar
pankmer index -t 10 --rounds 10 -g pankmer_input_fastas/only_HTF -o only_HTF.tar
```

## 1.2. Create adjacency matrices
```bash
pankmer adj-matrix -i ALL_genes.tar -o ALL_genes.adjmatrix.tsv
pankmer adj-matrix -i without_HTF.tar -o without_HTF.adjmatrix.tsv
pankmer adj-matrix -i only_HTF.tar -o only_HTF.adjmatrix.tsv

```

## 1.3. Create Heatmap with Jacard distances
```bash
pankmer clustermap -i ALL_genes.adjmatrix.tsv -o ALL_genes.adjmatrix.jaccard.svg --metric jaccard --width 20 --height 20
```

# 2. Create PCoA from adjacency matrices
I used the adjacency matrices from each of the datasets to compute Jaccard distances. I then used the distances to do PCoA's using the `pcoa` function from the `skikit-bio` python package.  

Additionally, I used the metainfomation to color the points in the PCoA figures.  

The analyses can be found in the Jupyter Notebook.

# 3. Create NJ trees
Using the Jaccard distances, I created NJ trees with the `DistanceTreeConstructor` function of the `Biophython` package. The code can be found at the same notebook and the produced trees are: `NJ_trees/ALL_genes.distance.jaccard.NJ.newick`, `NJ_trees/only_HTF.distance.jaccard.NJ.newick` and `NJ_trees/without_HTF.distance.jaccard.NJ.newick`

# 4. visualize NJ trees

i used the py script and the R to visual the PCoA and final NJ trees in two layouts, radial and polar.

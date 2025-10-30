#/bin/bash -l

#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=1:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/2025keykmer_74HTF
#$ -N keykmer
#$ -V
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/logs
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/logs
source ~/miniconda3/bin/activate phylogeny_snp
data=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs
mkdir $data/HTFrefs_with_wholegenome_7refs_rm_selfHTF
cp $data/HTF_7refs/*fa $data/HTFrefs_with_wholegenome_7refs_rm_selfHTF/
cp $data/excludeselfHTF/step3_wg.fa.excludeselfHTF/cleaned_noHTF/*fasta $data/HTFrefs_with_wholegenome_7refs_rm_selfHTF/
#then rename the .fasta to .fa


data2=$data/../forkmer_step1_kmerref_afterhamming

KMER=31
JFSIZE=300000000
wd=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.2_HTF_bykmers

cd $wd
mkdir -p intermediate_kmersfiltering intermediate_kmersfiltering/kmers_all intermediate_kmersfiltering/kmers_unique logs kmers_unique_hamming2

# --- Step 2: Count unique kmers per reference ---
echo "[2] Extracting unique kmers per reference"
for ref in ${data}/HTFrefs_with_wholegenome_7refs_rm_selfHTF/*.fa; do
    name=$(basename "$ref" .fa)

    jellyfish count -m $KMER -s 1000000 -t 10 -C -o ${wd}/intermediate_kmersfiltering/kmers_all/${name}.jf "$ref"
    jellyfish dump -c intermediate_kmersfiltering/kmers_all/${name}.jf | awk '$2==1{print $1}' > ${wd}/intermediate_kmersfiltering/kmers_all/${name}_uniquekmers.txt
done

# --- Step 3: Filter reference-specific kmers ---
echo "[3] Filtering allele-specific kmers"
for focal in ${wd}/intermediate_kmersfiltering/kmers_all/*_uniquekmers.txt; do
    name=$(basename "$focal" _uniquekmers.txt)
    
    # Build "others" file by concatenating all other kmer lists except $focal
    tmp_others="${wd}/tmp_others.txt"
    > "$tmp_others"  # clear file
    for other in ${wd}/intermediate_kmersfiltering/kmers_all/*_uniquekmers.txt; do
        [[ "$other" == "$focal" ]] && continue
        cat "$other" >> "$tmp_others"
    done

    # Sort and remove overlapping lines
    sort "$focal" > ${wd}/tmp_self.txt
    sort "$tmp_others" > tmp_others.sorted.txt
    
    comm -23 ${wd}/tmp_self.txt tmp_others.sorted.txt > ${wd}/intermediate_kmersfiltering/kmers_unique/${name}_unique.txt
done

# Clean up
rm ${wd}/tmp_self.txt ${wd}/tmp_others.txt ${wd}/tmp_others.sorted.txt


#then perform hamming...


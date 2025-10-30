#get the seq of each gene in the annotated p25.C2 from Sergio's annotation folder with provided coordinate in the gff annotation file from paper
#bash step1_getseq.sh &&
#use the seqs of genes query the p25.c2 renaed.contigs.reference and have the coordinate in chrom start end
#bash step2_minimaptop25c2.sh
#calculate the cov and dep (dep of genes/avg deg of whole ref) for each genes in 46 historical (40 ATUE5 + 6 nonATUE5) and 85 (55+30non) modern samples
#bash step3_cov_dep.sh
#bash step4_cov_dep_modern.sh
#summary the relative depth and plot a heatmap
bash step5_generate_table_relativedepth.sh
bash step6_generate_table_covprop.sh
#run binary matrix
bash step7_binary_table_75percent_on_both_depth_and_covprop.sh

#method2 see if we can assemble genes from mapped reads for these 46+95 reads
 

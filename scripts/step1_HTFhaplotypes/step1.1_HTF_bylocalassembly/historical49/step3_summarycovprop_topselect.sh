#!/bin/bash

# Directories
samples_txt="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/mappings/sampleswithassemblies.txt"
ls -d /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/assemblies/*/contigs.fasta 2>/dev/null | xargs -n1 dirname | xargs -n1 basename > $samples_txt
paf_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/results/step1_HTFhaplotypes/step1.1_HTF_bylocalassembly/historical49paf/mappings"
out_txt="$paf_dir/top_matched_haplotypes_by_propcoverage_detailed.txt"

# Output header
echo -e "sample\ttop_HTFhaplotype\tHTF_prop\ttfa_bestHTF_prop\ttop_TFA\tTFA_prop\tHTF_len\tHTF_cov\ttfa_bestHTF_len\ttfa_bestHTF_cov\tbestTFA_len\tbestTFA_cov" > "$out_txt"
# Loop through samples
while read -r sample; do
    paf="$paf_dir/${sample}_mapped.paf"
    [[ ! -f "$paf" ]] && continue

awk -v sname="$sample" '
{
    ref = $1
    qlen = $2
    qstart = $3
    qend = $4
    covlen = qend - qstart

    if (ref ~ /^HTF_/) {
        hap = ref
        gsub(/^HTF_/, "", hap)
        htf_cov[hap] += covlen
        if (qlen > htf_len[hap]) htf_len[hap] = qlen
    } else if (ref ~ /^TFA_/) {
        hap = ref
        gsub(/^TFA_/, "", hap)
        tfa_cov[hap] += covlen
        if (qlen > tfa_len[hap]) tfa_len[hap] = qlen
    }
}
END {
    best_htf = ""; best_tfa = ""; max_htf = 0; max_tfa = 0

    for (r in htf_cov) {
        if (htf_len[r] > 0) {
            prop = htf_cov[r] / htf_len[r]
            if (prop > max_htf) {
                max_htf = prop
                best_htf = r
            }
        }
    }

    for (r in tfa_cov) {
        if (tfa_len[r] > 0) {
            prop = tfa_cov[r] / tfa_len[r]
            if (prop > max_tfa) {
                max_tfa = prop
                best_tfa = r
            }
        }
    }

    if (best_htf != "") {
        htfp = htf_cov[best_htf] / htf_len[best_htf]
        hlen = htf_len[best_htf]
        hcov = htf_cov[best_htf]

        if (best_tfa != "") {
            tfap = tfa_cov[best_tfa] / tfa_len[best_tfa]
            tlen = tfa_len[best_tfa]
            tcov = tfa_cov[best_tfa]
        } else {
            tfap = 0; tlen = 0; tcov = 0
        }

        # TFA value corresponding to best HTF haplotype
        if (best_htf in tfa_cov) {
            tfabest_prop = tfa_cov[best_htf] / tfa_len[best_htf]
            tfabest_len = tfa_len[best_htf]
            tfabest_cov = tfa_cov[best_htf]
        } else {
            tfabest_prop = 0; tfabest_len = 0; tfabest_cov = 0
        }

        print sname "\t" best_htf "\t" htfp "\t" tfabest_prop "\t" best_tfa "\t" tfap "\t" hlen "\t" hcov "\t" tfabest_len "\t" tfabest_cov "\t" tlen "\t" tcov 
    }
}' "$paf" >> "$out_txt"

done < "$samples_txt"

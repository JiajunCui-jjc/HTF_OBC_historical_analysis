#wd = /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/2025_summer_paperfig_m57/data/HTFreference/forkmer_step0_7refs/excludeselfHTF/step2_HTFseqs_samestrandasinwg/

#p13.C7
#HTF_p5.D5	1803	12	1793	+	5	401475	18162	19943	1781	1781	60	tp:A:P	cm:i:181	s1:i:1781	s2:i:0	dv:f:0.0003	rl:i:0
echo 'p13.C7'
samtools faidx ../p13.C7.fasta "5:18151-19953" 
#double check if the two match
cat ../../refs/HTF_p13.C7.fa
#then write to sHTF* in $wd

#HTF_p21.F9	1245	4	1235	-	1	773043	752816	754047	1231	1231	60	tp:A:P	cm:i:117	s1:i:1231	s2:i:0	dv:f:0.0004	rl:i:0
echo 'p21.F9'
samtools faidx ../p21.F9.fasta "1:752807-754051" | seqtk seq -r -
cat ../../refs/HTF_p21.F9.fa

#HTF_p25.A12	1383	0	1373	-	7	24085	5838	7211	1373	1373	60	tp:A:P	cm:i:131	s1:i:1373	s2:i:0	dv:f:0.0004	rl:i:0
echo 'p25.A12'
samtools faidx ../p25.A12.fasta "7:5829-7211" | seqtk seq -r -

cat ../../refs/HTF_p25.A12.fa

#HTF_p25.C2	1803	12	1793	+	22	174356	18252	20033	1781	1781	60	tp:A:P	cm:i:175	s1:i:1781	s2:i:0	dv:f:0.0003	rl:i:0
echo 'p25.C2'
samtools faidx ../p25.C2.fasta "22:18241-20043" 

cat ../../refs/HTF_p25.C2.fa
#HTF_p26.D6	1803	12	1793	-	524	80855	9618	11399	1781	1781	60	tp:A:P	cm:i:184	s1:i:1781	s2:i:0	dv:f:0.0003	rl:i:0


echo 'p26.D6'
samtools faidx ../p26.D6.fasta "524:9609-11411"| seqtk seq -r - 

cat ../../refs/HTF_p26.D6.fa

#HTF_p23.B8	1803	12	1793	-	1	407267	387433	389214	1781	1781	60	tp:A:P	cm:i:181	s1:i:1781	s2:i:0	dv:f:0.0003	rl:i:0
echo 'p26.E7'
samtools faidx ../p26.E7.fasta "1:387424-389226"| seqtk seq -r -

cat ../../refs/HTF_p26.E7.fa

#HTF_p7.G11	1830	10	1820	-	59	234583	215654	217464	1810	1810	60	tp:A:P	cm:i:178	s1:i:1810	s2:i:0	dv:f:0.0003	rl:i:0
echo 'p7.G11'
samtools faidx ../p7.G11.fasta "59:215645-217474"| seqtk seq -r -

cat ../../refs/HTF_p7.G11.fa

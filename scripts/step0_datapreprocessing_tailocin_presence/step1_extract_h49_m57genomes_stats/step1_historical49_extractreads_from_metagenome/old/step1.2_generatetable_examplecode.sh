#!/bin/bash -l

#$ -pe mpi 5
#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=2:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
samplename=${1}

cd  /SAN/ugi/plant_genom/jiajucui/answers_h40/
all=$(cat answers_for_${samplename}.txt | grep 'q1 What' -A 5 | grep 'in total (Q' | cut -d '+' -f 1)
At=$(cat answers_for_${samplename}.txt | grep 'q1 What' -A 5 | grep 'mapped (1' | cut -d '+' -f 1)
res=`echo "scale=3; $At/$all" | bc`

allAt=$(cat answers_for_${samplename}.txt | grep 'q1.1 What' -A 4 | grep '119667750')
Atatleastoneread=$(cat answers_for_${samplename}.txt | grep 'base mapped to At ref:' -A 1 | cut -d ':' -f 2)
Atdepth=$(cat answers_for_${samplename}.txt | grep 'read depth to At ref:' -A 1 | sed -n 2p)

At_cov=`echo "scale=5; $Atatleastoneread/$allAt" | bc`

At_depth=`echo "scale=5; $Atdepth" | bc`

Ps=$(cat answers_for_${samplename}.txt | grep 'q2 What' -A 2 | grep 'mapped (1' | cut -d '+' -f 1)
unmapped=`echo "scale=3; $all-$At" | bc`
res2=`echo "scale=5; $Ps/$unmapped" | bc`

nonannotated=$(($unmapped-$Ps))
res3=`echo "scale=5; $nonannotated/$all" | bc`

allPs=$(cat answers_for_${samplename}.txt | grep 'q4 What' -A 4 | grep '5941411')
atleastoneread=$(cat answers_for_${samplename}.txt | grep 'bases with at least 1 read:' -A 1 | cut -d ':' -f 2)
depth=$(cat answers_for_${samplename}.txt | grep 'q6' -A 1 | sed -n 2p)

res4=`echo "scale=5; $atleastoneread/$allPs" | bc`

res5=`echo "scale=5; $depth" | bc`

echo $samplename 0$res 0$At_cov 0$At_depth 0$res2 0$res4 0$res5

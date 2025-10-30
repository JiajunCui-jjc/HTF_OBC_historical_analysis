echo 'samplename At_percent At_covered At_depth Ps_percent Ps_covered Ps_depth'  >> 2024_162_4th_answers_forplot.txt 


for line in $(cat samplenames.txt | sed -n $i'p'| awk '{print $1}');
do
bash plotdatagenerate.sh $line >> summarydata.txt;
done


To plot the length distribution of the 46 historical genomes (plants and pathogens):
We do check first the length distribution and found there is a peak at 75bp in pseudomonas genomes from Lopez 2022 (n=10).
So supp fig1 I include 29 PL and 7 HB samples that are clear in length distribution in both species. 
And include the Lopez 10 samples in supp fig2 where I calculate the insert size after this strategy of capturing reads:
First run adapter removal to collapse r1 and r2 (captured the molecular covered with overlapped reads).
And align r1.truncated and r2.truncated to the reference separately then merge them into one bam to capture reads that not overlapped but mapped to reference from r1 or r2, we suppose these reads are mostly 75bp length and if we calculate the molecular size by bam mapping info successfully, we can recover the length of molecular that were covered by two 75bp reads at two ends but not overlapped since the molecular is too long (longer than 150bp).

Then I combine them and calculate the length of these molecular as the length distribution of these genomes mapped to either plant or pathogen of the 10 Lopez samples.


The 27 and 109 seem like has a short avg length in both species but when I go back to the mapdamage length distribution there are a lot 75bp in pseudomonas genome of these two samples (you can see there are some reads around 150-250...)
So the 10 Lopez all contains big molecular in pseudomonas genomes.
While the 143 Lopez together give a beautiful negative regression between collapsed reads prop and insert size, the 10 together is neutral...

We could say that we saw all 10 samples have collapsed prop below around 65% so to maximum the info we have, for these reads we do an alternative strategy...
Or include the 46h correlation? But we dont have r1 and r2 for HB samples
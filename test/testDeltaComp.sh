# testDeltaComp.sh: Script for comparing the \delta
# computation by mlRho, Delta24, and the "true" values
# obtained directly from the template using the program
# delta. This is a faster and - I hope - more transparent
# test than that based on makedata.py.
# Date: August 7, 2014
# Author: Bernhard Haubold
LENGTH=100000
echo "Reference length: " $LENGTH
echo COMMENT Simulate diploid genome
ms 2 1 -p 10 -t $(( LENGTH / 100 )) -r $(( LENGTH / 200 )) $LENGTH | ms2dna > template.fasta
echo COMMENT Simulate genome sequencing
sequencer -a -c 10 -E 0.001 -P template.fasta > reads.fasta
getSeq -s mate=1 reads.fasta > mate1.fasta
getSeq -s mate=2 reads.fasta > mate2.fasta
echo COMMENT Run delta
delta -d 1 -D 5 template.fasta | cut -f 1,6
echo COMMENT Extract haploid template
getSeq -s S1 template.fasta > t.fasta
mv t.fasta template.fasta
echo COMMENT Map reads to haploid template using bowtie2
bowtie2-build --quiet template.fasta template 
bowtie2 -x template -p 4 -X 750 -f -S true-align.sam -1 mate1.fasta -2 mate2.fasta

#remove broken records
grep -v '\*' true-align.sam > true-align-clean.sam
mv true-align-clean.sam true-align.sam

samtools view -S -b true-align.sam 2> /dev/null > test.bam
echo COMMENT Prepare mlRho run
samtools sort test.bam test2 2> /dev/null
samtools index test2.bam 2> /dev/null
samtools mpileup 2> /dev/null test2.bam |
cut -f 2,5 | 
awk -f bam2pro.awk | 
formatPro > /dev/null
echo COMMENT Run mlRho
echo Distance   Delta
mlRho -M 0 -m 1 -M 5  | grep '^[1-9]' | awk '{print $1 "\t" $4}' | perl -pe 's/</\t/g' | awk '{print $1 "\t" $3}'
echo COMMENT Run Delta24
echo Distance   Delta
../src/Delta24 test2.bam 1 6 1 | grep '^D=' | perl -pe 's/=/\t/g;s/,//g' | awk '{print $2 "\t" $4}'

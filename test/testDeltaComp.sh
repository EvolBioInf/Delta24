# testDeltaComp.sh: Script for comparing the \delta
# computation by mlRho, Delta24, and the "true" values
# obtained directly from the template using the program
# delta. This is a faster and - I hope - more transparent
# test than that based on makedata.py.
# Date: August 7, 2014
# Author: Bernhard Haubold
#
# NB: when running bowtie, we need to map to the diploid template
# consisting of sequences S1 and S2. If we'd map to only one
# these two sequences, the genetic diversity of the resulting mapping
# would be halved. I am still puzzled why this is.
# But once we do map to the diploid template, we get
# TWO mappings, one to S1 and the other to S2. This doubles the
# length of our alignment. I merge the two alignments by
# (i) removing the @SQ header for S2
# (ii) replacing the template label S2 by S1
LENGTH=100000
echo "Reference length: " $LENGTH
echo COMMENT Simulate diploid genome
ms 2 1 -t $(( LENGTH / 10 )) -r $(( LENGTH / 200 )) $LENGTH | ms2dna > template.fasta
echo COMMENT Simulate genome sequencing
sequencer -a -c 10 -E 0.001 -P template.fasta > reads.fasta
getSeq -s mate=1 reads.fasta | sed -r 's/ f.+$//' > mate1.fasta
getSeq -s mate=2 reads.fasta | sed -r 's/ f.+$//' > mate2.fasta
echo COMMENT Map reads using bowtie
bowtie-build --quiet template.fasta template 
bowtie -p 4 -X 750 -f -S template -1 mate1.fasta -2 mate2.fasta 2> /dev/null > true-align.sam
# Remove @SQ header for S2
grep -v SN:S2 true-align.sam > t.sam
mv t.sam true-align.sam
# Replace the template label S2 by S1
#perl -pe 's/S2/S1/;s/ f.+$//' true-align.sam > t.sam
perl -pe 's/S2/S1/' true-align.sam > t.sam
mv t.sam true-align.sam
echo COMMENT Convert sequencer reads to BAM file
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
mlRho -m 1 -M 5  | grep '^[1-9]' | awk '{print $1 "\t" $4}' | perl -pe 's/</\t/g' | awk '{print $1 "\t" $3}'
echo COMMENT Run Delta24
echo Distance   Delta
../src/Delta24 test2.bam 1 6 1 | grep '^D=' | perl -pe 's/=/\t/g;s/,//g' | awk '{print $2 "\t" $4}'
echo COMMENT Run delta
delta -d 1 -D 5 template.fasta | cut -f 1,6

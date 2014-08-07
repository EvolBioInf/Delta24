echo COMMENT Simulate diploid genome
ms 2 1 -t 10000 -r 8000 1600000 | ms2dna > template.fasta
echo COMMENT Simulate genome sequencing
python makedata.py template.fasta 10 > /dev/null
echo COMMENT Convert sequencer reads to BAM file
samtools view -S -b true-align.sam 2> /dev/null > test.bam
echo COMMENT Prepare mlRho run
samtools sort test.bam test2 2> /dev/null
samtools index test2.bam 2> /dev/null
samtools view -b test2.bam contig_1:1-1000000 2> /dev/null | 
samtools mpileup - 2> /dev/null | 
cut -f 2,5 | 
awk -f bam2pro.awk | 
formatPro > /dev/null
echo COMMENT Run mlRho
echo Distance   Delta
mlRho -m 1 -M 5 | grep '^[1-9]' | awk '{print $1 "\t" $4}' | perl -pe 's/</\t/g' | awk '{print $1 "\t" $3}'
echo COMMENT Run Delta24
echo Distance   Delta
/usr/bin/time -f "%e %M" ../src/Delta24 test.bam 1 6 1 | grep '^D=' | perl -pe 's/=/\t/g;s/,//g' | awk '{print $2 "\t" $4}'
delta template.fasta -d 1 -D 5 | cut -f 1,6

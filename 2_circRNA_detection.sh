#!/usr/bin/env bash
# identify circRNAs from RNAseq 

echo ">>> aligning reads"
bowtie2 -p8 --very-sensitive --mm -M20 --score-min=C,-15,0 -x bt2_Homo_sapiens.GRCh38.dna.chromosome -q -U /media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/fastq/ASD_SRR5938456_1.fastq,/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/fastq/ASD_SRR5938456_2.fastq 2> circASD_SRR5938456_pass1.log | samtools view -hbuS - | samtools sort - circASD_SRR5938456
echo ">>> get the unmapped"
samtools view -hf 4 circASD_SRR5938456.bam | samtools view -Sb - > unmapped_circASD_SRR5938456.bam
echo ">>> split into anchors"
../unmapped2anchors.py unmapped_circASD_SRR5938456.bam | gzip > circASD_SRR5938456_anchors.fastq.gz
echo ">>> run find_circ.py"
mkdir circASD_SRR5938456_output
bowtie2 --reorder --mm -M20 --score-min=C,-15,0 -q -x bt2_Homo_sapiens.GRCh38.dna.chromosome -U circASD_SRR5938456_anchors.fastq.gz 2> circASD_SRR5938456_pass2.log | ../find_circ.py -G Homo_sapiens.GRCh38.dna.chromosome.1.fa -p circASD_SRR5938456_ -s circASD_SRR5938456_output/sites.log > circASD_SRR5938456_output/sites.bed 2> circASD_SRR5938456_output/sites.reads
echo ">>> Done."


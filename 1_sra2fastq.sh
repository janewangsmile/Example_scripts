#!/bin/bash
# convert file from sra to fastq

DATA=/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/sra/
FastQ=/media/jun/circRNAs/ASD_circRNAs_data/Brain_ASD_seq/fastq/

for i in $(ls $DATA)
do
	name=${i/'.sra'/''}
	~/Downloads/sratoolkit.2.10.0-ubuntu64/bin/fastq-dump -I --split-files $name -O $FastQ
done




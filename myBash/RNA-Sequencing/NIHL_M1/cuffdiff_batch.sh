#! /bin/bash

#cuffdiff -o No_Noise -b genome.fa -p 8 -L Apical,Basal -u merged_asm/merged.gtf \
#./tophat2_output/Barcode3/accepted_hits.bam,./tophat2_output/Barcode7/accepted_hits.bam \
#./tophat2_output/Barcode4/accepted_hits.bam,./tophat2_output/Barcode8/accepted_hits.bam 

#cuffdiff -o 90_dB -b genome.fa -p 8 -L Apical,Basal -u merged_asm/merged.gtf \
#./tophat2_output/Barcode1/accepted_hits.bam,./tophat2_output/Barcode13/accepted_hits.bam \
#./tophat2_output/Barcode2/accepted_hits.bam,./tophat2_output/Barcode14/accepted_hits.bam 

#cuffdiff -o 96_dB -b genome.fa -p 8 -L Apical,Basal -u merged_asm/merged.gtf \
#./tophat2_output/Barcode5/accepted_hits.bam,./tophat2_output/Barcode15/accepted_hits.bam \
#./tophat2_output/Barcode6/accepted_hits.bam,./tophat2_output/Barcode16/accepted_hits.bam 

#cuffdiff -o 96_dB_ZO -b genome.fa -p 8 -L Apical,Basal -u merged_asm/merged.gtf \
#./tophat2_output/Barcode9/accepted_hits.bam,./tophat2_output/Barcode11/accepted_hits.bam \
#./tophat2_output/Barcode10/accepted_hits.bam,./tophat2_output/Barcode12/accepted_hits.bam 

cuffdiff -o Basal_96_90 -b genome.fa -p 8 -L 96dB,90dB -u merged_asm/merged.gtf \
./tophat2_output/Barcode6/accepted_hits.bam,./tophat2_output/Barcode16/accepted_hits.bam \
./tophat2_output/Barcode2/accepted_hits.bam,./tophat2_output/Barcode14/accepted_hits.bam 

cuffdiff -o Basal_96ZO_90 -b genome.fa -p 8 -L 96dB_ZO,96dB -u merged_asm/merged.gtf \
./tophat2_output/Barcode10/accepted_hits.bam,./tophat2_output/Barcode12/accepted_hits.bam \
./tophat2_output/Barcode6/accepted_hits.bam,./tophat2_output/Barcode16/accepted_hits.bam 

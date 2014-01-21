#! /bin/bash


cuffdiff -o Basal_NoNoise_90dB -b NIHL_M1/genome.fa -p 8 -L NoNoise,90dB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode4/accepted_hits.bam,NIHL-5/tophat2_output/Barcode14/accepted_hits.bam,NIHL-6/tophat2_output/Barcode6/accepted_hits.bam,NIHL-6/tophat2_output/Barcode14/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode2/accepted_hits.bam,NIHL-5/tophat2_output/Barcode8/accepted_hits.bam,NIHL-6/tophat2_output/Barcode10/accepted_hits.bam,NIHL-6/tophat2_output/Barcode16/accepted_hits.bam 

cuffdiff -o Basal_NoNoise_96dB -b NIHL_M1/genome.fa -p 8 -L NoNoise,96dB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode4/accepted_hits.bam,NIHL-5/tophat2_output/Barcode14/accepted_hits.bam,NIHL-6/tophat2_output/Barcode6/accepted_hits.bam,NIHL-6/tophat2_output/Barcode14/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode6/accepted_hits.bam,NIHL-5/tophat2_output/Barcode16/accepted_hits.bam,NIHL-6/tophat2_output/Barcode2/accepted_hits.bam,NIHL-6/tophat2_output/Barcode8/accepted_hits.bam 

cuffdiff -o Basal_NoNoise_96dBZO -b NIHL_M1/genome.fa -p 8 -L NoNoise,96dBZO -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode4/accepted_hits.bam,NIHL-5/tophat2_output/Barcode14/accepted_hits.bam,NIHL-6/tophat2_output/Barcode6/accepted_hits.bam,NIHL-6/tophat2_output/Barcode14/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode10/accepted_hits.bam,NIHL-5/tophat2_output/Barcode12/accepted_hits.bam,NIHL-6/tophat2_output/Barcode4/accepted_hits.bam,NIHL-6/tophat2_output/Barcode12/accepted_hits.bam 

cuffdiff -o Basal_90dB_96dB -b NIHL_M1/genome.fa -p 8 -L 90dB,96dB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode2/accepted_hits.bam,NIHL-5/tophat2_output/Barcode8/accepted_hits.bam,NIHL-6/tophat2_output/Barcode10/accepted_hits.bam,NIHL-6/tophat2_output/Barcode16/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode6/accepted_hits.bam,NIHL-5/tophat2_output/Barcode16/accepted_hits.bam,NIHL-6/tophat2_output/Barcode2/accepted_hits.bam,NIHL-6/tophat2_output/Barcode8/accepted_hits.bam 

cuffdiff -o Basal_90dB_96dBZO -b NIHL_M1/genome.fa -p 8 -L 90dB,96dBZO -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode2/accepted_hits.bam,NIHL-5/tophat2_output/Barcode8/accepted_hits.bam,NIHL-6/tophat2_output/Barcode10/accepted_hits.bam,NIHL-6/tophat2_output/Barcode16/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode10/accepted_hits.bam,NIHL-5/tophat2_output/Barcode12/accepted_hits.bam,NIHL-6/tophat2_output/Barcode4/accepted_hits.bam,NIHL-6/tophat2_output/Barcode12/accepted_hits.bam 

cuffdiff -o Basal_96dB_96dBZO -b NIHL_M1/genome.fa -p 8 -L 96dB,96dBZO -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode6/accepted_hits.bam,NIHL-5/tophat2_output/Barcode16/accepted_hits.bam,NIHL-6/tophat2_output/Barcode2/accepted_hits.bam,NIHL-6/tophat2_output/Barcode8/accepted_hits.bam \
NIHL-5/tophat2_output/Barcode10/accepted_hits.bam,NIHL-5/tophat2_output/Barcode12/accepted_hits.bam,NIHL-6/tophat2_output/Barcode4/accepted_hits.bam,NIHL-6/tophat2_output/Barcode12/accepted_hits.bam 


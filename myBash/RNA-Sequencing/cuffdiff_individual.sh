#! /bin/bash

cuffdiff -o NoNoise_20251_2LB_vs_NE90_20251_1LB  -b NIHL_M1/genome.fa  -p 8 -L NoNoise_20251_2LB,NE90_20251_1LB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode4/accepted_hits.bam NIHL-5/tophat2_output/Barcode2/accepted_hits.bam 

cuffdiff -o NE96dB_20251_3LB_vs_NE90_20252_2LB  -b NIHL_M1/genome.fa  -p 8 -L NE96dB_20251_3LB,NE90_20252_2LB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode6/accepted_hits.bam NIHL-5/tophat2_output/Barcode8/accepted_hits.bam 

cuffdiff -o NE96dBZO_20254_1LB_vs_NE96ZO_20255_5LB  -b NIHL_M1/genome.fa  -p 8 -L NE96dBZO_20254_1LB,NE96ZO_20255_5LB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode10/accepted_hits.bam NIHL-5/tophat2_output/Barcode12/accepted_hits.bam 

cuffdiff -o NoNoise_20256_2LB_vs_NE96_20257_3LB  -b NIHL_M1/genome.fa  -p 8 -L NoNoise_20256_2LB,NE96_20257_3LB -u merged_asm/merged.gtf \
NIHL-5/tophat2_output/Barcode14/accepted_hits.bam NIHL-5/tophat2_output/Barcode16/accepted_hits.bam 

cuffdiff -o NE96dB_20241_2LB_vs_NE96dBZO_20241_3LB  -b NIHL_M1/genome.fa  -p 8 -L NE96dB_20241_2LB,NE96dBZO_20241_3LB -u merged_asm/merged.gtf \
NIHL-6/tophat2_output/Barcode2/accepted_hits.bam NIHL-6/tophat2_output/Barcode4/accepted_hits.bam 

cuffdiff -o NoNoise_20245_2RB_vs_NE96dB_20245_3LB  -b NIHL_M1/genome.fa  -p 8 -L NoNoise_20245_2RB,NE96dB_20245_3LB -u merged_asm/merged.gtf \
NIHL-6/tophat2_output/Barcode6/accepted_hits.bam NIHL-6/tophat2_output/Barcode8/accepted_hits.bam 

cuffdiff -o NE90dB_20245_5RB_vs_NE96dBZO_20251_4LB  -b NIHL_M1/genome.fa  -p 8 -L NE90dB_20245_5RB,NE96dBZO_20251_4LB -u merged_asm/merged.gtf \
NIHL-6/tophat2_output/Barcode10/accepted_hits.bam NIHL-6/tophat2_output/Barcode12/accepted_hits.bam 

cuffdiff -o NoNoise_20257_1RB_vs_NE90dB_20257_2LB  -b NIHL_M1/genome.fa  -p 8 -L NoNoise_20257_1RB,NE90dB_20257_2LB -u merged_asm/merged.gtf \
NIHL-6/tophat2_output/Barcode14/accepted_hits.bam NIHL-6/tophat2_output/Barcode16/accepted_hits.bam 



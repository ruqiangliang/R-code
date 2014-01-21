#!/bin/bash

paste -d '' <(echo; sed -n '1,${n;p;}' lane2_NoIndex_L002_R2_001.fastq | sed G) lane2_NoIndex_L002_R1_001.fastq | sed '/^$/d' | fastx_barcode_splitter.pl --bol --mismatches 1 --bcfile ../NIHL_M1/barcode.txt --prefix debarcoded --suffix ".fq"

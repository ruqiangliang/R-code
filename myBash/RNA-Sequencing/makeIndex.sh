#!/bin/bash

echo "This bash script make bam file index *.bmi"
for f in $(ls); do
	cd $f
	samtools index accepted_hits.bam
	cd ..
done

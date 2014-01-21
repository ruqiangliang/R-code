#!/bin/bash

for f in $(ls); do
	cd $f
	samtools index accepted_hits.bam
	cd ..
done

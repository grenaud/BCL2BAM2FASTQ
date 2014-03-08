all: 
	make -C libgab
	make -C bam2fastq
	make -C BCL2BAM
	make -C fastq2bam


clean:
	make -C libgab clean
	make -C bam2fastq clean
	make -C BCL2BAM clean
	make -C fastq2bam clean


.PHONY: all

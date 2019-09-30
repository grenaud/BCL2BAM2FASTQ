all: bamtools/build/src/api/libbamtools.a libgab/utils.o
	make -C libgab
	make -C bam2fastq
	make -C BCL2BAM
	make -C fastq2bam

libgab/utils.h:
	rm -rfv libgab/
	mkdir libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git libgab/


libgab/utils.o: bamtools/build/src/api/libbamtools.a  libgab/utils.h
	make -C libgab/

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	mkdir bamtools/
	git clone  https://github.com/pezmaster31/bamtools.git bamtools/ && bamtools/ &&  git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

bamtools/build/src/api/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


clean:
	make -C libgab clean
	make -C bam2fastq clean
	make -C BCL2BAM clean
	make -C fastq2bam clean


.PHONY: all

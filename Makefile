all: bamtools/build/src/api/libbamtools.a libgab/libgab.a
	make -C libgab
	make -C bam2fastq
	make -C BCL2BAM
	make -C fastq2bam

libgab/libgab.h:
	rm -rfv libgab/
	mkdir libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git libgab/


libgab/libgab.a: bamtools/build/src/api/libbamtools.a  libgab/libgab.h
	make -C libgab/


bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	mkdir bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git bamtools/ && cd bamtools/ #&& git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

bamtools/build/src/api/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ &&  mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi  && make && cd ../..


clean:
	make -C libgab clean
	make -C bam2fastq clean
	make -C BCL2BAM clean
	make -C fastq2bam clean


.PHONY: all

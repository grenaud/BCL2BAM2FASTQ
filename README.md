About
----------------------

BCL2BAM2FASTQ is a series of programs to convert from :

1) BCL to BAM

2) BAM to fastq

3) fastq to BAM

Downloading:
----------------------

Go to https://github.com/grenaud/BCL2BAM2FASTQ and either:

1) Download ZIP 

or

2) Do a "git clone --recursive https://github.com/grenaud/BCL2BAM2FASTQ.git"


Installation:

----------------------

1) Build Bamtools first:

cd bamtools/
mkdir build/
cd build/
cmake ..
make 
cd ../..

2) Build the submodules and main code by typing :

make



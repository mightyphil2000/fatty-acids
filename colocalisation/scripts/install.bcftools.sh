wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
bunzip2 bcftools-1.11.tar.bz2 
tar -xvf bcftools-1.11.tar 
cd bcftools-1.11
./configure --prefix=/mnt/storage/home/ph14916



cd samtools-1.x    # and similarly for bcftools and htslib
./configure --prefix=/where/to/install
make
make install

The executable programs will be installed to a bin subdirectory under your specified prefix, so you may wish to add this directory to your $PATH:

export PATH=/mnt/storage/home/ph14916/bin:$PATH    # for sh or bash users

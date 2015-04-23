Cheat sheet for Unix: [http://www.rain.org/~mkummel/unix.html](http://www.rain.org/~mkummel/unix.html)

Tutorials for Unix: [http://www.ee.surrey.ac.uk/Teaching/Unix/](http://www.ee.surrey.ac.uk/Teaching/Unix/)


###  To go on the cluster:
ssh jfortin@jhpce02.jhsph.edu

### TO create a directory:
mkdir <nameOfYourDirectory>

### To access directory:
cd <nameOfYourDirectory>

### To download a file (from GEO):
curl -O ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP014%2FSRP014134/SRR518875/SRR518875.sra


### To use sra-toolkit installed in my session:
/home/student/jfortin/sratoolkit.2.4.5-2-centos_linux64/bin/fastq-dump

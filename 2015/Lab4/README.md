Cheat sheet for Unix: [http://www.rain.org/~mkummel/unix.html](http://www.rain.org/~mkummel/unix.html)

Tutorials for Unix: [http://www.ee.surrey.ac.uk/Teaching/Unix/](http://www.ee.surrey.ac.uk/Teaching/Unix/)


To go on the cluster:

    ssh jfortin@jhpce02.jhsph.edu

Basic commands:

    mkdir <nameOfYourDirectory>
    cd <nameOfYourDirectory>

### Some downloads:

    curl -O ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP014%2FSRP014134/SRR518875/SRR518875.sra

To download SRA toolkit:

    curl -O http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.5-2/sratoolkit.2.4.5-2-ubuntu64.tar.gz

and unzip it:

    tar -zxvf sratoolkit.2.4.5-2-ubuntu64.tar.gz
 
Let's set up the permissions:

    chmod 777 sratoolkit.2.4.5-2-ubuntu64/bin/fastq-dump


Let's get on a node:

    qrsh -l mem_free=1G
    cd homework3

To download the fastq file:

    ../sratoolkit.2.4.5-2-ubuntu64/bin/fastq-dump SRR518875.sra

To get the number of reads:

    wc -l SRR518875.fastq (Divided by 4)


We will create a subset of the data: (1000 first reads):

    head -n 4000 SRR518875.fastq > example.fastq



Let's download the Yeast genome:

    curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/s_cerevisiae.ebwt.zip
    unzip s_cerevisiae.ebwt.zip

### Alignment using bowtie
    module load bowtie
    bowtie ../yeast_genome/s_cerevisiae example.fastq example.sam 
    
# Suppose you wanna get rid of columns 5,6,7,8:
bowtie  --suppress 5,6,7,8 ../yeast_genome/s_cerevisiae example.fastq example.sam 






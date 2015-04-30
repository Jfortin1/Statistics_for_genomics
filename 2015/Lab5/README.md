# Lab 5: Analysis of ChIP-Seq data

## Part 1: Using macs

Macs is a simple and powerfull model-based peak calling algorithm for ChIP-Seq data on the command line. On the cluster, it is already installed and you can load it by

    module load macs
    
If you want to have more information about how macs works, here is a nice tutorial:

[http://www.slideshare.net/lucacozzuto/macs-course](http://www.slideshare.net/lucacozzuto/macs-course)


To go on the cluster:

    ssh jfortin@jhpce02.jhsph.edu

Basic commands:

    mkdir <nameOfYourDirectory> # To create a directory 
    cd <nameOfYourDirectory> # To access a directory
    touch <newFile> # To create a new file

### Some downloads:

To download the SRA experiment in the homework3 folder:
    
    cd homework3
    curl -O ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP014%2FSRP014134/SRR518875/SRR518875.sra

Also need to download SRA toolkit that will allow us to convert the SRA experiment to a fastq file:

    cd $home
    curl -O http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.5-2/sratoolkit.2.4.5-2-ubuntu64.tar.gz
    tar -zxvf sratoolkit.2.4.5-2-ubuntu64.tar.gz
 
Set up permissions for fastq-dump:

    chmod 777 sratoolkit.2.4.5-2-ubuntu64/bin/fastq-dump

### Converting the SRA experiment to fastq file:

We first get on a node using the qrsh command:

    qrsh -l mem_free=1G
    
We are ready to extract the fastq file from the SRA experiment:
    cd homework3

    cd homework3
    ../sratoolkit.2.4.5-2-ubuntu64/bin/fastq-dump SRR518875.sra

We now have the file SRR518875.fastq in the homework directory. 

    head -n 8 SRR518875.fastq 
    wc -l SRR518875.fastq # Number of reads x 4


For this lab we will work only with a subset of the data (first 1000 reads): 

    head -n 4000 SRR518875.fastq > example.fastq


### Alignment using bowtie

We need the Yeast genome index files from the Bowtie website:

    curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/s_cerevisiae.ebwt.zip
    unzip s_cerevisiae.ebwt.zip

We are ready to align the reads with the default parameters of Bowtie:
    
    module load bowtie
    bowtie ../s_cerevisiae example.fastq example.sam 
    
To get rid of columns 5,6,6,8 in the output:

    bowtie  --suppress 5,6,7,8 ../yeast_genome/s_cerevisiae example.fastq example.sam 






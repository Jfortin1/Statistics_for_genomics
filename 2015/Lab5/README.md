# Lab 5: Analysis of ChIP-Seq data

## Part 1: Using macs

Macs is a simple and powerfull model-based peak calling algorithm for ChIP-Seq data on the command line. On the cluster, it is already installed and you can load it by

    module load macs
    
If you want to have more information about how macs works, here is a nice tutorial:

[http://www.slideshare.net/lucacozzuto/macs-course](http://www.slideshare.net/lucacozzuto/macs-course)

### Example dataset for this lab:

I'm going to use the first million reads for the treatment and control samples (first condition) from the GSE39147 experiment. I stored the reads in the files `treatment.fastq` and `control.fastq`.    

First step is to align the reads using `bowtie`. The option `--sam` will output sam files:

    module load bowtie
    bowtie --sam s_cerevisiae/s_cerevisiae treatment.fastq treatment.sam 
    bowtie --sam s_cerevisiae/s_cerevisiae control.fastq control.sam 
    
Make sure you are poitning to the right genome path (see Lab 4). Let's use macs to get the peaks:

    module load macs
    macs14 -t treatment.sam -c control.sam -g 1.26e+7 Condition1

`-t` specifies the input file. macs accepts different formats: sam, bam, and bed for instance. The `-c` specficies the control samples. The `-g` specifies the length of the genome. How do I know that s. cerevisiae has length 1.26e+7? I used google :). THe output files will have names Condition1_(andSomethingElse)




    






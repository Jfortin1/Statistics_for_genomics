# Lab 5: Analysis of ChIP-Seq data

## Part 1: Using macs to produce the list of peaks

Macs is a simple and powerfull model-based peak calling algorithm for ChIP-Seq data on the command line. On the cluster, it is already installed and you can load it by

    module load macs
    
If you want to have more information about how macs works, here is a nice tutorial:

[http://www.slideshare.net/lucacozzuto/macs-course](http://www.slideshare.net/lucacozzuto/macs-course)

### Example dataset for this lab:

I'm going to use the first million reads for the treatment and control samples (first condition) from the GSE39147 experiment. Note that the experiment has two conditions and you are asked to compare the peaks between the two conditions. I stored the first million reads for the first condtion in the files `treatment.fastq` and `control.fastq`.    

First step is to align the reads using `bowtie`. The option `--sam` will output sam files:

    module load bowtie
    bowtie --sam s_cerevisiae/s_cerevisiae treatment.fastq treatment.sam 
    bowtie --sam s_cerevisiae/s_cerevisiae control.fastq control.sam 
    
Make sure you are poitning to the right genome path (see Lab 4). Let's use macs to get the peaks:

    module load macs
    macs14 -t treatment.sam -c control.sam -g 1.26e+7 -n Condition1 -w

`-t` specifies the input file. macs accepts different formats: sam, bam, and bed for instance. The `-c` specficies the control samples. The `-g` specifies the length of the genome. How do I know that s. cerevisiae has length 1.26e+7? I used google :). `-n` specifices the output files; the output files will have names `Condition1_(andSomethingElse)`. Finally, the option `-w` will create wiggle files for each chromosome analyzed. 

The ouput files are `Condition1_negative_peaks.xls`, `Condition1_peaks.bed`, `Condition1_peaks.xls` and `Condition1_summits.bed`.

The file `peaks.bed` contains the peaks information:

    head Condition1_peaks.bed
    
    Scchr02 44810   45956   MACS_peak_1     124.31
    Scchr02 89067   90214   MACS_peak_2     142.75
    Scchr02 167630  168780  MACS_peak_3     173.45
    Scchr02 332053  333140  MACS_peak_4     260.70
    Scchr02 414863  416113  MACS_peak_5     255.51
    Scchr02 477800  479149  MACS_peak_6     60.38
    Scchr02 592433  593497  MACS_peak_7     131.19
    Scchr02 603820  604934  MACS_peak_8     190.47
    Scchr02 605509  606660  MACS_peak_9     209.67
    Scchr03 10      783     MACS_peak_10    77.97

The columns are: chr, start, end, peak id and score (-10*log(pvalue)). The file `summits.bed` contains the summits of the peaks:

    head Condition1_summits.bed
    
    Scchr02 45393   45394   MACS_peak_1     59.00
    Scchr02 89648   89649   MACS_peak_2     62.00
    Scchr02 168194  168195  MACS_peak_3     76.00
    Scchr02 332617  332618  MACS_peak_4     110.00
    Scchr02 415576  415577  MACS_peak_5     93.00
    Scchr02 478442  478443  MACS_peak_6     50.00
    Scchr02 593003  593004  MACS_peak_7     84.00
    Scchr02 604346  604347  MACS_peak_8     92.00
    Scchr02 606035  606036  MACS_peak_9     97.00
    Scchr03 203     204     MACS_peak_10    62.00
    
The columns are: chr, start, end, peak id and score (height of the summit). 
    
## Part 2: Annotation of the peaks

[ChIPpeakAnno](http://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) and [ChIPseeker](http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) are two useful packages for annotation of ChIP-Seq peaks. They are available via Bioconductor.

    source("http://bioconductor.org/biocLite.R")
    biocLite("ChIPpeakAnno")
    biocLite("ChIPseeker")

I won't go through all the details for these two packages -- the vignettes are very well documented and you should read them for your homework. Here are some tips to get started:
    
    # To read the peaks into R:
    peaks <- read.csv("Condition1_peaks.xls", sep="\t", skip=23)
    head(peaks)
    
              chr  start    end length summit tags X.10.log10.pvalue. fold_enrichment
    1 Scchr02  44811  45956   1146    584  154             124.31            4.14
    2 Scchr02  89068  90214   1147    582  156             142.75            4.52
    3 Scchr02 167631 168780   1150    565  166             173.45            5.54
    4 Scchr02 332054 333140   1087    565  185             260.70            7.91
    5 Scchr02 414864 416113   1250    714  200             255.51            6.77
    6 Scchr02 477801 479149   1349    643  183              60.38            2.67
      FDR...
    1   4.29
    2   4.69
    3   6.12
    4  14.29
    5  13.33
    6   2.63
    
We will create a GenomicRanges object that contain the genomic intervals of the peaks and all the other columns above. One way to do that is to use the `data.frame2GRanges` function from the `bsseq` package:

    data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
        stopifnot(class(df) == "data.frame")
        stopifnot(all(c("start", "end") %in% names(df)))
        stopifnot(any(c("chr", "seqnames") %in% names(df)))
        if("seqnames" %in% names(df))
            names(df)[names(df) == "seqnames"] <- "chr"
        if(!ignoreStrand && "strand" %in% names(df)) {
            if(is.numeric(df$strand)) {
                strand <- ifelse(df$strand == 1, "+", "*")
                strand[df$strand == -1] <- "-"
                df$strand <- strand
            }
            gr <- GRanges(seqnames = df$chr,
                          ranges = IRanges(start = df$start, end = df$end),
                          strand = df$strand)
        } else {
            gr <- GRanges(seqnames = df$chr,
                          ranges = IRanges(start = df$start, end = df$end))
        }
        if(keepColumns) {
            dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                         "DataFrame")
            elementMetadata(gr) <- dt
        }
        names(gr) <- rownames(df)
        gr
    }

    library(GenomicRanges)
    peaks.gr <- data.frame2GRanges(peaks)
    elementMetadata(peaks.gr) <- peaks[,4:9]
    
To use the ChIPseeker package, you will need the annotation `TxDb.Scerevisiae.UCSC.sacCer3.sgdGene` instead of the annotation `TxDb.Hsapiens.UCSC.hg19.knownGene` used in the vignette examples. 

    source("http://bioconductor.org/biocLite.R")
    biocLite("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
    library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
    
As an example, here is how to produce a coverage plot of the peaks using the log10 pvalue as a score:
    
    library(ChIPseeker)
    covplot(peaks.gr, weightCol="X.10.log10.pvalue.")

## Part 3: Visualization of the peaks usint the wig files

    library(rtracklayer)
    wigDir <- "Condition1_MACS_wiggle/treat"
    track  <- import(file.path(wigDir, "Condition1_treat_afterfiting_Scchr01.wig.gz"))
    pos  <- (start(track)+end(track))/2
    score <- track$score
    indices <- 1000:1500
    barplot( score[indices], col="firebrick")

## Part 4: More about the Hopkins Cluster

- ssh, qstat, qmem, qdel, screen -S, qdel -u, qsub -cwd -V -l mem_free=10G,h_vmem=12G yourScript.sh $i;-l h_fsize=40G

To set up ssh keys:

[Click here](https://jhpce.jhu.edu/knowledge-base/authentication/login/)


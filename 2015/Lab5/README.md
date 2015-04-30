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

[ChIPpeakAnno](http://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) and `ChIPseeker` are two useful packages for annotation of ChIP-Seq peaks. They are available via Bioconductor [http://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html](http://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) and [http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html]()http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html

    source("http://bioconductor.org/biocLite.R")
    biocLite("ChIPpeakAnno")

I won't go through all the details of the package -- the vignette is very well documented and I suggest reading it for your homework: [http://www.bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.pdf](http://www.bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.pdf). Here's how to get started:

    library(ChIPpeakAnno)
    macsOutput <- read.csv("Condition1_negative_peaks.xls", sep="\t")
    head(macsOutput)
    
       chr    start   end length summit tags X.10.log10.pvalue. fold_enrichment
    1 Scmito     7  3251   3245    776  941            1284.94            5.40
    2 Scmito 65578 68862   3285   1896  762             800.00            6.07
    3 Scmito 68893 71886   2994    386  547             209.30            4.84
    
Let's create a GenomicRanges object. We first need this useful function from the `bsseq` package:

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



    library(bsseq)
    macsGR <- bsseq::data.frame2GRanges(macsOutput)
    elementMetadata(macsGR) <- macsOutput[,4:8]
    
## Part 3: More about the Hopkins Cluster



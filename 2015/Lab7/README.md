# Lab 7: RNA-Seq analysis

We will follow a lab from HarvardX PH525x Data Analysis for Genomics that can be found  [Here](https://github.com/genomicsclass/labs/blob/master/course5/rnaseq_gene_level.Rmd)

## Part 1: 

[Click here](link)


```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("airway")
library(airway)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
```

```{r}
library(airway)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
csv.file <- file.path(dir, "sample_table.csv")
sample.table <- read.csv(csv.file, row.names=1)
bam.files <- file.path(dir, paste0(sample.table$Run, "_subset.bam"))
gtf.file <- file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf")
```

```{r}
```


```{r}
library(Rsamtools)
bam.list <- BamFileList(bam.files)
library(GenomicFeatures)
# for Bioc 3.0 use the commented out line
# txdb <- makeTranscriptDbFromGFF(gtf.file, format="gtf")
txdb <- makeTxDbFromGFF(gtf.file, format="gtf")
exons.by.gene <- exonsBy(txdb, by="gene")
```






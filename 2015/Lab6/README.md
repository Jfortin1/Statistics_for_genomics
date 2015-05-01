# Lab 6

## Part 1: Feedback on previous assignments

- For each plot, make sure to label your axes correctly, have a meaningful title, and describe what the plot is supposed to tell us in the caption. 
- For trimming reads with Bowtie , one can use the following options as part of the `bowtie` command:
```
--trim3 10 --trim5 10 
```
which will trims the reads by 10 bps respectively for the 3' and 5' ends. 
- To produce a SAM file, use the option `--sam'
- To produce a fastq file containing the unaligned reads, add
```
--un alignedOutput
```
which will save the unaligned reads in the file alignedOutput.fastq.
- While a SAM file is a regular text file, the BAM file is a binary version of the SAM file. One can convert a SAM file to a BAM files by using `samtools` already installed on the cluster:
```
module load samtools
samtools view -b myfile.sam -o myfile.bam
```
- For homework 3, you were asked to compute the coverage for each chromosome. The coverate can be calculated by counting how many reads you have for each chromosome, multiply by the length of the reads (36) and divide by the length of the chromosome. There are many ways to do that by either using the command line or R. First, let me introduce some useful commands in Unix.
- `cut` is a command line utility very useful to access columns in a text file. For instance, consider the following file `example.txt` containing three columns and 7 lines: 
```
chr1	123123	123124
chr2	12323	123123
chr1	143123	153124
chr2	12	24
chr2	12	24
chr3	133	1212
chr4	999	1024
```
If you want to grab the second column, you can do:
```
cut -f 2 example.txt

123123
12323
143123
12
12
133
999
```
Now, if you want to compute how many times each chromsome appears, you can do
```
cut -f 1 example.txt | sort | uniq -c

   2 chr1
   3 chr2
   1 chr3
   1 chr4
```
First, the `|` operator (usually refered as "pipe") will use the output of the left-hand side command as the input of the right-hand side command. In the example above, `cut -f 1 example.txt` will grad the first column, than it will be sorted using `sort`, and the sorted list will be fed to `uniq -c`. The latter command will create a list a unique elements, and the option `-c` will output the  number of occurences for each element. That's one way to count how many reads mapped to each of the chromosome, for instance. 


## Part 2: More about the Hopkins Cluster

- ssh, qstat, qmem, qdel, screen -S, qdel -u, qsub -cwd -V -l mem_free=10G,h_vmem=12G yourScript.sh $i;-l h_fsize=40G

To set up ssh keys:

[Click here](https://jhpce.jhu.edu/knowledge-base/authentication/login/)




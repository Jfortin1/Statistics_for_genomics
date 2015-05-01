# Lab 6

## Part 1: Feedback on previous assignments

- For each plot, make sure to label your axes correctly, have a meaningful title, and describe what the plot is supposed to tell us in the caption. 
- For homework 3, you were asked to compute the coverage for each chromosome. The coverate can be calculated by counting how many reads you have for each chromosome, multiply by the length of the reads (36) and divide by the length of the chromosome. 
- For trimming reads with Bowtie , one can use the following options as part of the `bowtie` command:
```
--trim3 10 --trim5 10 
```
which will trims the reads by 10 bps respectively for the 3' and 5' ends. 
  





## Part 2: More about the Hopkins Cluster

- ssh, qstat, qmem, qdel, screen -S, qdel -u, qsub -cwd -V -l mem_free=10G,h_vmem=12G yourScript.sh $i;-l h_fsize=40G

To set up ssh keys:

[Click here](https://jhpce.jhu.edu/knowledge-base/authentication/login/)




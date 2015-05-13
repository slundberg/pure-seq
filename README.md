# pure-seq
A program that learns how to remove background signal from ChIP-seq data using hundreds of control tracks.

One idea for the eventual usage is: (just to get some guess out there of what it will be)

```bash
pure-seq input.bam > output.sam
```

Each read in the output would be labeled with the background poisson arrival rate at that position (say the starting position of the read). This could then be used by downstream processing algorithms to weight each read appropreitely. In fact we could also provide a simple peak caller that just thresholded windows of reads by their poisson p-value:

```bash
poisson-peak output.sam --pvalue 0.001 --window 1000 > peaks.bed
```


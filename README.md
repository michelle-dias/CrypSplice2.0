# CrypSplice 2.0

## Overview
CrypSplice 2.0 is a robust computational tool designed to detect and analyze cryptic alternative splicing events from RNA-seq data. It provides enhanced capabilities for identifying differential splicing patterns that are not annotated (cryptic) in standard gene models.

## Installation
git clone https://github.com/michelle-dias/CrypSplice2.0 
pip install .

## Usage
### Cryptic Junctions 
crypsplice CrypticJunctions \
  -c1 "./chr19_Bams/Contrl1.sorted.bam,./chr19_Bams/Contrl2.sorted.bam,./chr19_Bams/Contrl3.sorted.bam" \
  -c2 "./chr19_Bams/Treat1.sorted.bam,./chr19_Bams/Treat2.sorted.bam,./chr19_Bams/Treat3.sorted.bam" \
  -gtf /path/to/annotations.gtf \
  -fasta /path/to/genome.fa \
  -s 2 \
  -annotated \
  -o ./CrypticJunctions/

### Cryptic Load - Diff
crypsplice CrypticLoad Diff \
  -c1 "./chr19_Bams/Contrl1.sorted.bam,./chr19_Bams/Contrl2.sorted.bam,./chr19_Bams/Contrl3.sorted.bam" \
  -c2 "./chr19_Bams/Treat1.sorted.bam,./chr19_Bams/Treat2.sorted.bam,./chr19_Bams/Treat3.sorted.bam" \
  -gtf $gtf \
  -fasta $fasta \
  -s 2 \
  -o ./CrypticLoad_Diff/


### Cryptic Load - Clust
crypsplice CrypticLoad Clust \
  -samples "./chr19_Bams/Contrl1.sorted.bam,./chr19_Bams/Contrl2.sorted.bam,./chr19_Bams/Contrl3.sorted.bam,./chr19_Bams/Treat1.sorted.bam,./chr19_Bams/Treat2.sorted.bam,.chr19_Bams/Treat3.sorted.bam \
  -p $processors \
  -gtf $gtf \
  -fasta $fasta \
  -s 2 \
  -o ./CrypticLoad_Clust/ 


### Cryptic Load - Clust
crypsplice CrypticLoad Clust \
  -samples "./chr19_Bams/Contrl1.sorted.bam,./chr19_Bams/Contrl2.sorted.bam,./chr19_Bams/Contrl3.sorted.bam,./chr19_Bams/Treat1.sorted.bam,./chr19_Bams/Treat2.sorted.bam,.chr19_Bams/Treat3.sorted.bam \
  -p $processors \
  -gtf $gtf \
  -fasta $fasta \
  -s 2 \
  -o ./CrypticLoad_Clust/ 


### PlotJunctions
crypsplice PlotJunctions \
    -pj Novel_Junctions_output \
    -c1 ./chr19_Bams/Contrl1.sorted.bam,./chr19_Bams/Contrl2.sorted.bam,./chr19_Bams/Contrl3.sorted.bam \
    -c2 ./chr19_Bams/Treat1.sorted.bam,./chr19_Bams/Treat2.sorted.bam,.chr19_Bams/Treat3.sorted.bam \
    -gtf $gtf \
    -fasta $fasta \ 
    -s 2 \
    -o ./PlotJunctions/


## Arguments / Options
CrypticJunctions – Infer cryptic splice junctions across two conditions
    Required arguments:
          -c1 : Control BAM files (space-separated list)
          -c2 : Treated BAM files (space-separated list)
          -gtf : Reference GTF annotation
          -fasta : Reference FASTA file
          -s : Strand specificity (0 = unstranded, 1 = frd, 2 = rev)
          -o : Output directory
    Optional arguments:
          -prefix : Output file prefix (default: CrypSplice)
          -p : Number of processors to use (default: 10)
          -annotated : Include annotated junctions
          -j : Read count cutoff (default: 10)
          -l : Junction length cutoff (default: 50)
          -pa_p : pOverA filter: P value (default: 0.5)
          -pa_a : pOverA filter: A value (default: 15)
          -g : Gene filter (remove multigene junctions)
          -b : Batch file (TSV: Sample\tBatch)

CrypticLoad – Calculate and analyze cryptic load

      Subcommand: Diff
          Required arguments:
                -c1 : Control BAM files (space-separated list)
                -c2 : Treated BAM files (space-separated list)
                -gtf : Reference GTF annotation
                -fasta : Reference FASTA file
                -s : Strand specificity (0 = unstranded, 1 = frd, 2 = rev)
                -o : Output directory
          Optional arguments:
                -prefix : Output file prefix (default: CrypSplice)
                -p : Number of processors (default: 10)
                -b : Batch file (TSV: Sample\tBatch)
                -j : Junction read count cutoff (default: 10)
                -d : Original count cutoff (default: 0)
                -l : Intron length cutoff (default: 50)
                -pa_p : pOverA filter: P value (default: 0.5)
                -pa_a : pOverA filter: A value (default: 0)

      Subcommand: Clust
          Required arguments:
              -samples : List of cryptic load files (space-separated)
              -gtf : Reference GTF annotation
              -fasta : Reference FASTA file
              -s : Strand specificity (0 = unstranded, 1 = frd, 2 = rev)
              -o : Output directory
          Optional arguments:
              -r : Manually specify cluster rank
              -i : NMF iterations to infer rank (default: 50)
              -prefix : Output file prefix (default: CrypSplice)
              -p : Number of processors (default: 10)
              -b : Batch file (TSV: Sample\tBatch)
              -j : Junction read count cutoff (default: 10)
              -d : Original count cutoff (default: 0)
              -l : Intron length cutoff (default: 50)
              -pa_p : pOverA filter: P value (default: 0.5)
              -pa_a : pOverA filter: A value (default: 0)


PlotJunctions – Visualize cryptic junctions of interest 
      Required arguments:
          -o : Output directory
          -pj : Junctions file (Novel_Junctions.txt or Annotated_Junctions.txt)
          -c1 : Control BAM files (space-separated list)
          -c2 : Treated BAM files (space-separated list)
          -gtf : Reference GTF annotation
          -fasta : Reference FASTA file
          -s : Strand specificity (0 = unstranded, 1 = frd, 2 = rev)
          -b : Batch file (optional)
      Optional arguments:
          -prefix : Output file prefix (default: CrypSplice)
          -p : Number of processors (default: 10)
          -top_n : Number of top junctions to plot (default: 10)
          -j_ids : Specific junction IDs to plot (comma-separated)
          --merge : Merge bigWigs for each condition (C1 and C2)
          -c1_names : Names for C1 replicates or group name if merging
          -c2_names : Names for C2 replicates or group name if merging
          -c1_color : Hex color code for C1 tracks (default: #2387de)
          -c2_color : Hex color code for C2 tracks (default: #e04196)
          -j : Read count cutoff (default: 10)
          -l : Junction length cutoff (default: 50)


## Requirements
Python >=3.9.5
numpy >=1.26.4
sys >=3.9.5
pandas >=1.5.0
concurrent.futures >=3.9.5
pybedtools >=0.8.2
subprocess >=3.9.5
sklearn >=1.1.2
matplotlib >=3.6.21
gc >=3.9.5
time >=3.9.5
pyBigWig >=0.3.18




## License
This project is licensed under the MIT License – see the ./LICENSE file for details 


## Contact
Author: Michelle Dias
Email: michelledias10@gmail.com
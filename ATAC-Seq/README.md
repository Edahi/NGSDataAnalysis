# ATAC-Seq
Assay for Transposase-Accessible Chromatin followed by Sequencing focus on the obtention of genomic regions that are not being wrapped by either nucleosomes or other cellular elements, meaning that they are "Accessible" for interaction with cellular elements.

### 0- (Optional) Download SRR fastq data files
- Download the compressed FASTQ files using [SRA-toolkit's](https://github.com/ncbi/sra-tools) fastq-dump.

### 1- Quality check
- Sequencing lanes are merged before quality assessment
- General sequencing quality, length, barcodes and different oddities are first analyzed through [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
<img src="ATAC-Seq/Images/00.QualityControl.png" width="200">
On the left figure is the data you DON'T want to ideally work with. It is totally worth discussing with your peers about what could have happened that caused the high noise.
The left figure illustrates the importance of visualizing your raw data.
This particularly important if you are working with public data. It is likely that the authors submitted the raw data and therefore did not remove the barcodes as the one illustrated here.
Dececting these variations early helps to transition smoothly to the Differential Accessibility analysis.

### 2- Step mapping
_1st-step Mapping_
- _Raw_ reads are aligned to reference genome (e.g. mm10) using [Botwie1](http://bowtie-bio.sourceforge.net/manual.shtml)
- Unmapped reads are processed with [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) through the use of [trim_galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Low-quality bases are removed from 5'/3' and as minimum, 1bp is clipped from the 3'-end.

_2nd-step Mapping_
- _Processed_ reads are mapped to the reference genome with more permissive parametrs, such a higher insertion distance and slight more mismatches rate.
- For comparison, _unmapped_ reads are quality-checked with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Both mapping results are merged and sorted using [samtools](http://samtools.sourceforge.net/).
<img src="ATAC-Seq/Images/01.2StepMapping.png" width="300" height="100">
Ocassionally, it might be questionable to use this strategy. Personally I think it is nor painstaking neither really time consuming and the benefits widely varies depending the samples themselves.
Illustrated here is the **Gain** of performing the 2-Step mapping, on the first row, there is close to 5% gain, on the second row, I gained close to 17% (the most common scenario) and on the third one I gained almost double the information, 90%.
The reason to align the data w/o any initial trimming or filtering is because not all reads needs it. Let map all the raw data and apply the filter to those that didn't make it.

### 3- Cleaning
- Reads aligning to the mitochodrial chromosome are eliminated through the use of both [samtools](http://samtools.sourceforge.net/) and a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script.
- Reads mapping to [ENCODE's blackliested regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) are removed with [bedtools](https://bedtools.readthedocs.io/en/latest/).
- PCR duplicates are removed with [Picard-tools](https://broadinstitute.github.io/picard/)

There are two strategies to perform ATAC-Seq.
1. Using a buffer to reduce Mitochondrial contamination and 
2. Don't.
which are [omniATAC-Seq]() and [ATAC-Seq]() respectively. The advantages of using one versus the other are better illustrated in my figure below. 
<img src="ATAC-Seq/Images/04.Filtering.png">
Each row illustrates the mapping results and then the reads that reamins after the respective filters.
The final two columns indicate the data loss and the methodology used to obtain the sequencing reads.
The top three rows had a loss of **above 50% mapped reads**. The bottom four had a loss of **less than 8%**. The difference was the methodology.

### 4- Fragment-length estimation
- Filtered reads are used to calculate the distance among read pairs. You may wish to contact [Xi Chen](mailto:xi.chen.xchen@gmail.com)
<img src="ATAC-Seq/Images/05.FragmentLengthEstimation.png">
Distance between ATAC-Seq mapped reads is expected to illustrate a pattern of nucleosome occupancy ( ≤ 200 bp wraps around nucleosomes). 
Should you have data that looks more like a flat slope, you might need to discuss biological insights of this ibservation to discard unproper data preparation.

### 5- Subcategories obtention
- Subnucleosomal reads are obtained with a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script (reads closer than 100bp).
- Tn5 insertion sites are obtained through the use of a perl script

### 6- Peak calling
- [HOMER](http://homer.ucsd.edu/homer/ngs/tagDir.html) tag directories are generated for the filtered reads and the subnucleosomal pairs. 
- Using these tagDirectories, peaks are called using the [findPeaks](http://homer.ucsd.edu/homer/ngs/peaks.html) parameters " -style dnase -region -nfr ".


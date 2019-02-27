# ATAC-Seq
Assay for Transposase-Accessible Chromatin followed by Sequencing focus on the obtention of genomic regions that are not being wrapped by either nucleosomes or other cellular elements, meaning that they are "Accessible" for interaction with cellular elements.

### 0- (Optional) Download SRR fastq data files
- Download the FASTQ files using [SRA-toolkit's](https://github.com/ncbi/sra-tools) fastq-dump.

### 1- Quality check
- Sequencing lanes are merged before quality assessment
- General sequencing quality, length, barcodes and different oddities are first analyzed through [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
<img src="/ATAC-Seq/Images/00.QualityControl.png">
On the left figure is the data you DON'T want to ideally work with. It is totally worth discussing with your peers about what could have happened that caused the high noise.
The middle graph is thehow average well-done data usually looks like. 
The left figure illustrates the importance of visualizing your raw data, which is particularly important if you are working with data of public origins ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi), etc.
It is likely that the authors submitted the raw data and therefore did not remove the barcodes as the one illustrated here.
Early detection of these variations usually leads to a smooth transition to the Differential Accessibility analysis.

### 2- Step mapping
_1st-step Mapping_
- _Raw_ reads are aligned to reference genome (e.g. mm10) using [Botwie1](http://bowtie-bio.sourceforge.net/manual.shtml)
- Unmapped reads are processed with [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) through the use of [trim_galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). 
Low-quality bases are removed from 5'/3' and as minimum, 1bp is clipped from the 3'-end.

_2nd-step Mapping_
- _Processed_ reads are mapped to the reference genome with more permissive parametrs, such a higher insertion distance and slight more mismatches rate.
- For comparison, _unmapped_ reads are quality-checked with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Both mapping results are merged and sorted using [samtools](http://samtools.sourceforge.net/).
<img src="/ATAC-Seq/Images/01.2StepMapping.png">
It might be questionable to use this strategy all of the time. Personally I think it is nor painstaking neither time consuming and the benefits widely varies depending the samples themselves.
Illustrated here is the *Gain* of performing the 2-Step mapping, on the first row, there is close to 5% gain, on the second row, I gained close to 17% (the most common scenario) and on the third one I gained almost double the information (90%) than that of mapping the raw reads only.
The reason to align the data without any initial trimming or filtering is because not all reads needs it. Let map all the raw data and apply the filter to those that didn't make it.

### 3- Contamination assessment
- Unmapped reads are assembled with [velvet and velveth](https://www.ebi.ac.uk/~zerbino/velvet/) into contigs.
- The longest and more covered _Top 3_ contigs are [_Blast_-ed](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to assess any kind of contamination on the sample.
<img src="/ATAC-Seq/Images/02.ContaminationAnalysis.png">
You can definitely skip this if your mapping results are around and above the 70%. Is your data is way below that, it might be worthy to further investigate possible conntamination sources.
On the middle row the sample is oddily unmapping to the reference genome. I took this one and proceed as described to determine that this sample was effectively contaminated with salmon DNA.
<img src="/ATAC-Seq/Images/03.ContaminationAnalysisIllustration.png">
It is important to discuss results with your peers so you can narrow down all of the multiple possibilities.

### 4- Cleaning
- Reads aligning to the mitochodrial chromosome are eliminated through the use of both [samtools](http://samtools.sourceforge.net/) and a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script.
- Reads mapping to [ENCODE's blackliested regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) are removed with [bedtools](https://bedtools.readthedocs.io/en/latest/).
- PCR duplicates are removed with [Picard-tools](https://broadinstitute.github.io/picard/)

There are two strategies to perform ATAC-Seq.
1. Using a buffer to reduce Mitochondrial material and 
2. Don't.

which are [omniATAC-Seq](https://www.nature.com/articles/nmeth.4396) and [ATAC-Seq](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142727.mb2129s109) respectively. 
The advantages of using one versus the other are better illustrated in my figure below. 
<img src="/ATAC-Seq/Images/04.Filtering.png">
Each row illustrates the mapping results and then the reads that reamins after the respective filters.
The final two columns indicate the data loss and the methodology used to obtain the sequencing reads.
The top three rows had a loss of __above 50% mapped reads__. The bottom four had a loss of __less than 8%__. The difference was the methodology.

### 5- Fragment-length estimation
- Filtered reads are used to calculate the distance among read pairs. You may wish to contact [Xi Chen](mailto:xi.chen.xchen@gmail.com)
<img src="/ATAC-Seq/Images/05.FragmentLengthEstimation.png">
Distance between ATAC-Seq mapped reads is expected to illustrate a pattern of nucleosome occupancy ( ≤ 200 bp wraps around nucleosomes). 
Should you have data that looks more like a flat slope, you might need to discuss biological insights of this ibservation to discard unproper data preparation.

### 6- Subcategories obtention
- Subnucleosomal reads are obtained with a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script (reads closer than 100bp).
- Tn5 insertion sites are obtained through the use of a perl script

### 7- Generate Genome-Browser tracks
- From the tagDirectories [makeUCSCfile](http://homer.ucsd.edu/homer/ngs/ucsc.html) is used to auto-generate a bedGraph.gz file with a factor size of 1e20.
- Decompress bedGraph.gz, head-trim and sorted by chromose and position to generate a proper bedGraph file.
- UCSC's [bedGraphToBigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html)'d the proper bedGraph to generate a bigWig file.
- A proper track line describing the sample and its location on the server is generated.

### 8- Peak calling
- [HOMER](http://homer.ucsd.edu/homer/ngs/tagDir.html) tag directories are generated for the filtered reads and the subnucleosomal pairs. 
- Using these tagDirectories, peaks are called using the [findPeaks](http://homer.ucsd.edu/homer/ngs/peaks.html) parameters " -style dnase -region -nfr ".
- Accessible peaks from the subnucleosomal reads data are also obtained with [MACS2](https://github.com/taoliu/MACS) parameters "-g mm -q 0.0001  --keep-dup all --nomodel --call-summits".
- MACS2 peaks summits (deconvolve subpeaks within each peak) are extend 200 bp each and resulting peaks are merged among samples (Master set)
<img src="/ATAC-Seq/Images/06.DefineAccessibleRegions.png">
The longest white lines above the smaller are the initial peaks called by MACS2, by obtaining the smallest lines, called summit (re-analysing peak-shape and defining "sub"-peaks) we can perform analysis with greater granularity.

## Version changes:

### (4.0.0) -- 2019/02/19
- Added Blast analysis using no mappable reads (Contamination Assessment)
- Change the way I calculated the raw total reads
- Change remapping parameters to extend search space
- Added the calculation of Usable Reads
- Modified the way I gather the summary results
- Moved the Peak calling section to the end
- Changed the way I generate the bigWigs files (Now I used my custom script)
- Modified How I call peaks with MACS2 to concentrate on the sumimts
- Expand the summits to 200bp from the center
- Changed email inside script

### (3.0.0) -- 2018/10/22
- Corrected minor filename conveniences (Add HOMER prefix to HOMER peaks)
- Corrected missing backslash to make final file consistent with template format
- Added MACS2 program and peaks
- Renamed folder "06.HOMER_Peaks" by "06.Peaks" given MACS2 peaks introduction.

### (2.0.0) -- 2018/06/28
- Data downloaded now is compressed
- Added summary statistics to track mapping results
- Changed the order between removing duplicates and filtering blacklisted regions
- Renamed TagAlign to Fragment Length Estimate
- Added code necessary to generate BigWigs
- Added option to remove intermediate files at the end
- Add Genome Browser tracks || Experimental with Tn5_9bp
- Added code to check if SRR was given

# NGSDataAnalysis

Welcome to my NGS Data Analysis script folder!
During my PhD in Bioinformatics and Systems Biology at the University of California San Diego I was exposed to very diverse sequencing technologies routinarily. 
Therefore it was an ideal senario that allowed me to write code for the automatization of the data analysis of most of these sequencing datasets.
On this repository I have added the bash scripts that I wrote for each sequencing technology.

## Sequencing folders

- [ATAC-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/ATAC-Seq)
- [ChIP-like-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/ChIP-like-Seq)
- [CMS-IP-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/CMS-IP-Seq)
- [GRO-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/GRO-Seq)
- [HMCP-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/HMCP-Seq)
- [RNA-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/RNA-Seq)
- [WGB-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/WGB-Seq)
- [WG-Seq](https://github.com/Edahi/NGSDataAnalysis/tree/master/WG-Seq)

## Overview of each sequencing technique

ATAC-Seq
--------
_Assay for Transposase-Accessible Chromatin followed by Sequencing_ yield genomic regions that are not wrapped by either nucleosomes or cellular elements, meaning that they are accessible for interaction with other molecules.

<img src="/ATAC-Seq/Images/07.DARs.png">

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/ATAC-Seq)

### Analysis  Overview:

1. Generate Directories Framework
2. Merge lanes & rm barcodes (if any).
3. FASTQC for downloaded file
4. 1st mapping (raw reads)
5. _trim\_galore_ for unmapped reads
6. FASTQC for cleaned Unmapped reads
7. 2nd mapping (filtered unmapped reads)
8. Get unmapped reads from 2nd round.
9. De-novo assembly from unmapped
10. Retrieval of top 3 longest assemblies
11. Blast on top three longest assemblies
12. Merging of mapping results
13. Remove chrM
14. Remove Duplicates
15. Remove Encode's Blacklist
16. Gather all mapping summaries
17. Extract Sub-Nucleosomal fragments
18. Obtain Tn5 footprint
19. Generate HOMER's TagDirectory
20. Generate BigWigs for SubNuc, Tn5 and whole mapping results
21. Add BigWigs to tracks file
22. Estimate fragment size
23. Explore PhantomPeaks
24. Call Peaks fromo SubNuc
25. Extend to 200bp the peak's summit
26. Remove Intermediate Files 

RNA-Seq
-------
RNA-Seq allows a cuantitative approach for the analysis of transcripts in bulk RNA from cells.

<img src="/RNA-Seq/Images/07.MDS.png" width="700" height="300">

<img src="/RNA-Seq/Images/08.InteractiveDEGexploration.png" width="700" height="300">

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/RNA-Seq)

### Analysis  Overview:

___Pre-Differential Expression Analysis___:

1. Generate Directories Framework
2. Check if Paired-End data
3. Obtain read length (Determine STAR reference)
4. Merge lanes & rm barcodes (if any)
5. FASTQC for downloaded file
6. Map reads to tRNAs and rRNA (contamination estimation).
7. STAR mapping
8. Generate HOMER's TagDirectory
9. Run RSeQC's analyses
10. Generate Multi Wig Tracks
11. Subreads' Feature Counts.
12. Mapping stats recollection

___Differential Expression Analysis___:

1. Load multiple libraries and counts data.
2. Establish conditions (from provided Config file)
3. Obtain EntrezID for future annotation and add to main dataset
4. Data Normalization:
5. Visualization of the normalization proof and change effect.
6. Interactive Samples Clustering 
7. Differential Analysis Design
8. Voom-Normalization (with graphs)
9. Applyance of the model per pairs
10. Heatmaps of the most Variable Genes
11. Interactive Differential Expression Plots
12. Genome Onthology 
13. GO over-representation Upregulated


ChIP-like-Seq
-------------
_Chromatin Immunoprecipitation followed by Sequencing_ (ChIP-Seq) and similar technologies obtains genomic regions associated to a protein or modification of interest.

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/ChIP-like-Seq)

### Analysis  Overview:

1. Generate Directories Framework
2. check if Paired-End data
3. Merge lanes & rm barcodes (if any)
4. FASTQC for downloaded file
5. 1st mapping
6. Filter unmapped reads
7. _trim\_galore_ for unmapped reads
8. FASTQC for cleaned Unmapped reads
9. 2nd mapping
10. Get mapped reads
11. Merging of both mapping results
12. Remove chrM
13. Remove Duplicates
14. Remove Encode's Blacklist
15. Gather all mapping summaries
16. Generate HOMER's TagDirectory
17. Create BigWigs
18. Add BigWigs to tracks file:
19. Estimate fragment size
20. Explore PhantomPeaks
21. Call Peaks
22. Remove intermediary results

CMS-IP-Seq
----------
_Cytosine-5-MethyleneSulfonate followed by Immunoprecipitation and Sequencing_ is a method that allows the obtention of genomic regions enriched for the 5hydroxymethylcytosine epigenetic modification.

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/CMS-IP-Seq)

### Analysis  Overview:

1. Generate Directories Framework
2. check if Paired-End data
3. FASTQC for downloaded file
4. BSmap mapping
5. Separate lambda from genomic reads
6. For each, calculate Conversion Efficiency
7. Generate BAM files
8. Remove Duplicates
9. Remove Encode's Blacklist
10. Estimate fragment size
11. Explore PhantomPeaks
12. Get the methylation balls
13. MACS2 Call Peaks w/o input

GRO-Seq
-------
_Global Run On followed by Sequencing_ allows the analysis of nacent RNA, that is, genes that are being transcribed at a certain time point.

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/GRO-Seq)

### Analysis  Overview:

1. Generate Directories Framework
- Currently functional for Single-End data
2. Merge lanes & rm barcodes (if any)
3. Remove Poly-A tails
4. FASTQC for downloaded file
5. [BWA] Alignment 
6. Generate HOMER's TagDirectory
7. Mapping stats recollection:

The reason this particular pipeline is so simple is because after the HOMER tag directory, the downstream analysis really varies depending on the biological question.

HMCP-Seq
---------
_HydroxyMethylCytocine immunoPrecipitation followed by Sequencing_ This is a method that also allows the obtention of genomic regions enriched for the 5hydroxymethylcytosine epigenetic modification. 
The difference with [CMS-IP](https://github.com/Edahi/NGSDataAnalysis#cms-ip-seq) are:

- The use of Spike-Ins for accurate 5hmC measurement
- The immunoprecipitation is through biotin instead of the CMS adduct.
- Reads alphabet is not reduced because there is not bisulfite treatment
- Alignment to the unmodified reference genome can be done straightforward with e.g. [Bowtie(2)](http://bowtie-bio.sourceforge.net/manual.shtml) or [BWA](http://bio-bwa.sourceforge.net/).

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/HMCP-Seq)

### Analysis  Overview:

As described above, this data is then analyzed with the [ChIP-like pipeline](https://github.com/Edahi/NGSDataAnalysis#chip-like-seq) with the addition of Spike-Ins to the ref genome.
The steps below detail the normalization method particular to this technique:

1. Separate HMCP Spike-in BAMs
2. Separate T4 Spike-in BAM
3. Separate hg38 (or mm10) Genome BAM
4. Generate windowed BigWigs [MEDIPs] for T4 and reference genome
5. Spike-In Normalization Factor w/INPUT data
	1. For Input and Reference, calculate the RPM normalization factor and keep the reference's RPM factor apart.
	2. Convert Input's and reference's BigWig coverage to RPM.
	3. Subtract Input's from reference's RPM.
	4. Reverse RPM scaling using the stored reference's RPM factor
	5. Calculate the T4 Spike-in Exo-scaling factor
	6. Scale Reference with the calculated Exo-factor
6. Export the resulting BigWig for further analysis
7. Report Total counts per step (sanity check).
8. [MACS2] Call peaks in the Spiked-In normalized data.


WGB-Seq
-------

_Whole-Genome Bisulfite followed by Sequencing_ is a marvelous and expensive technique that allows the study of methylation analysis at base resolution. It cannot differentiate between 5hmC and 5mC since both almost no not interact with the Sodium Bisulphite treatment the DNA is subjected to.
Therefore, unless OxBS-Seq or TAB-Seq (none of them discussed here... yet) or another similar technique that allows the base-resolution identification of 5hmC (or 5mC) is used, we have to refer to the WGBS signal as __5mC + 5hmC__.


[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/WGB-Seq)

### Analysis  Overview:

1. Generate Directories Framework
- Currently functional for Paired-End data only (haven't seen a SE WGBS data yet).
2. Merge data files
3. FASTQC files
4. BSmap mapping
5. Separate lambda (and phiX if used) from genomic reads.
6. For each, calculate Bisulfite treatment efficiency
7. Generate BAM files
8. Get the methylation balls
-This program offer the option to remove duplicates and overhang fragments.
9. Filter CG-context methylation levels
10. Downstream analysis

WG-Seq
------
_Whole-Genome Sequencing_ allows for the analysis of single nucleotide variants, chromosomal rearrangements and much more.

[Extended Analysis and code](https://github.com/Edahi/NGSDataAnalysis/tree/master/WG-Seq)

### Analysis  Overview:

1. Generate Directories Framework
2. Merge lanes & rm barcodes (if any)
3. FASTQC for downloaded file
4. 1st [BWA] aligment (mapping).
5. Filter unmapped reads
6. _trim\_galore_ for unmapped reads
7. FASTQC for cleaned Unmapped reads
8. 2nd [BWA] aligment (mapping).
9. Merging of both mapping results
10. Remove Duplicates
11. Remove Encode's Blacklist
12. Remove Random chromosomes
13. Estimate fragment size
14. Explore PhantomPeaks
15. Generate HOMER's TagDirectory
16. Create BigWigs
17. Mapping stats recollection

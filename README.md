#### NGSDataAnalysis
Welcome to my NGS Data Analysis script folder!
During my PhD in Bioinformatics and Systems Biology at the University of California San Diego I was exposed to very diverse sequencing technologies rounarily. 
Therefore it was an ideal experience that allowed me to write code for the automatization of the data analysis of most of these sequencing datasets.
On this repository I have added the bash scripts for multiple of these NGS.

ATAC-Seq
--------
_Assay for Transposase-Accessible Chromatin followed by Sequencing_ yield genomic regions that are not wrapped by either nucleosomes or cellular elements, meaning that they are accessible for interaction with other molecules.

Analysis  Overview:

#### 1. Generate Directories Framework
#### 2. Merge lanes & rm barcodes (if any).
#### 3. FASTQC for downloaded file
#### 4. 1st mapping (raw reads)
#### 5. _trim\_galore_ for unmapped reads
#### 6. FASTQC for cleaned Unmapped reads
#### 7. 2nd mapping (filtered unmapped reads)
#### 8. Get unmapped reads from 2nd round.
#### 9. De-novo assembly from unmapped
#### 10. Retrieval of top 3 longest assemblies
#### 11. Blast on top three longest assemblies
#### 12. Merging of mapping results
#### 13. Remove chrM
#### 14. Remove Duplicates
#### 15. Remove Encode's Blacklist
#### 16. Gather all mapping summaries
#### 17. Extract Sub-Nucleosomal fragments
#### 18. Obtain Tn5 footprint
#### 19. Generate HOMER's TagDirectory
#### 20. Generate BigWigs for SubNuc, Tn5 and whole mapping results
#### 21. Add BigWigs to tracks file
#### 22. Estimate fragment size
#### 23. Explore PhantomPeaks
#### 24. Call Peaks fromo SubNuc
#### 25. Extend to 200bp the peak's summit
#### 26. Remove Intermediate Files 

ChIP-like-Seq
-------------
_Chromatin Immunoprecipitation followed by Sequencing_ (ChIP-Seq) and similar technologies obtains genomic regions associated to a protein or modification of interest.

Analysis overview:

#### 1. Generate Directories Framework
#### 2. check if Paired-End data
#### 3. Merge lanes & rm barcodes (if any)
#### 4. FASTQC for downloaded file
#### 5. 1st mapping
#### 6. Filter unmapped reads
#### 7. _trim\_galore_ for unmapped reads
#### 7. FASTQC for cleaned Unmapped reads
#### 8. 2nd mapping
#### 9. Get mapped reads
#### 10. Merging of both mapping results
#### 11. Remove chrM
#### 12. Remove Duplicates
#### 13. Remove Encode's Blacklist
#### 14. Gather all mapping summaries
#### 15. Generate HOMER's TagDirectory
#### 16. Create BigWigs
#### 17. Add BigWigs to tracks file:
#### 18. Estimate fragment size
#### 19. Explore PhantomPeaks
#### 20. Call Peaks
#### 21. Remove intermediary results

RNA-seq
-------
RNA-Seq allows a cuantitative approach for the analysis of transcripts in bulk RNA from cells.

Analysis overview:

_Pre-Differential Expression Analysis_.

#### 1. Generate Directories Framework
#### 2. Check if Paired-End data
#### 3. Obtain read length (Determine STAR reference)
#### 4. Merge lanes & rm barcodes (if any)
#### 5. FASTQC for downloaded file
#### 6. Map reads to tRNAs and rRNA (contamination estimation).
#### 7. STAR mapping
#### 8. Generate HOMER's TagDirectory
#### 9. Run RSeQC's analyses
#### 10. Generate Multi Wig Tracks
#### 11. Subreads' Feature Counts.
#### 12. Mapping stats recollection

_Differential Expression Analysis_.

#### 1. Load multiple libraries and counts data.
#### 2. Establish conditions (from provided Config file)
#### 3. Obtain EntrezID for future annotation and add to main dataset
#### 4. Data Normalization:
#### 5. Visualization of the normalization proof and change effect.
#### 6. Interactive Samples Clustering 
#### 7. Differential Analysis Design
#### 8. Voom-Normalization (with graphs)
#### 9. Applyance of the model per pairs
#### 10. Heatmaps of the most Variable Genes
#### 11. Interactive Differential Expression Plots
#### 12. Genome Onthology 
#### 13. GO over-representation Upregulated

CMS-IP-Seq
----------
_Cytosine-5-MethyleneSulfonate followed by Immunoprecipitation and Sequencing_ is a method that allows the obtention of genomic regions enriched for the 5hydroxymethylcytosine epigenetic modification.

Analysis overview:

#### 1. Generate Directories Framework
#### 2. check if Paired-End data
#### 3. FASTQC for downloaded file
#### 4. BSmap mapping
#### 5. Separate lambda from genomic reads
#### 6. For each, calculate Conversion Efficiency
#### 7. Generate BAM files
#### 8. Remove Duplicates
#### 9. Remove Encode's Blacklist
#### 10. Estimate fragment size
#### 11. Explore PhantomPeaks
#### 12. Get the methylation balls
#### 13. MACS2 Call Peaks w/o input

WG-Seq
------
_Whole-Genome Sequencing_ allows for the analysis of single nucleotide variants, chromosomal rearrangements and much more.

Analysis overview:

#### 1. Generate Directories Framework
#### 2. Merge lanes & rm barcodes (if any)
#### 3. FASTQC for downloaded file
#### 4. 1st [BWA] aligment (mapping).
#### 5. Filter unmapped reads
#### 6. _trim\_galore_ for unmapped reads
#### 7. FASTQC for cleaned Unmapped reads
#### 8. 2nd [BWA] aligment (mapping).
#### 9. Merging of both mapping results
#### 10. Remove Duplicates
#### 11. Remove Encode's Blacklist
#### 12. Remove Random chromosomes
#### 13. Estimate fragment size
#### 14. Explore PhantomPeaks
#### 15. Generate HOMER's TagDirectory
#### 16. Create BigWigs
#### 17. Mapping stats recollection


HMCP-Seq
---------
_HydroxyMethylCytocine immunoPrecipitation followed by Sequencing_ This is a method that also allows the obtention of genomic regions enriched for the 5hydroxymethylcytosine epigenetic modification. 
The difference with [CMS-IP](https://github.com/Edahi/NGSDataAnalysis#cms-ip-seq) are:

- The use of Spike-Ins for accurate 5hmC measurement
- The immunoprecipitation is through biotin instead of the CMS adduct.
- Reads alphabet is not reduced because there is not bisulfite treatment
- Alignment to the unmodified reference genome can be done straightforward with e.g. [Bowtie(2)](http://bowtie-bio.sourceforge.net/manual.shtml) or [BWA](http://bio-bwa.sourceforge.net/).

Analysis overview:

As described above, this data is then analyzed with the [ChIP-like pipeline](https://github.com/Edahi/NGSDataAnalysis#chip-like-seq) with the addition of Spike-Ins to the ref genome.
The steps below detail the normalization method particular to this technique:

#### 1. Separate HMCP Spike-in BAMs
#### 2. Separate T4 Spike-in BAM
#### 3. Separate hg38 (or mm10) Genome BAM
#### 4. Generate windowed BigWigs [MEDIPs] for T4 and reference genome
#### 5. Spike-In Normalization Factor w/INPUT data
- For Input and Reference, calculate the RPM normalization factor and keep the reference's RPM factor apart.
- Convert Input's and reference's BigWig coverage to RPM.
- Subtract Input's from reference's RPM.
- Reverse RPM scaling using the stored reference's RPM factor
- Calculate the T4 Spike-in Exo-scaling factor
- Scale Reference with the calculated Exo-factor

#### 6. Export the resulting BigWig for further analysis
#### 7. Report Total counts per step (sanity check).
#### 8. [MACS2] Call peaks in the Spiked-In normalized data.


GRO-Seq
-------
_Global Run On followed by Sequencing_ allows the analysis of nacent RNA, that is, genes that are being transcribed at a certain time point.

Analysis overview:

#### 1. Generate Directories Framework
- Currently functional for Single-End data

#### 2. Merge lanes & rm barcodes (if any)
#### 3. Remove Poly-A tails
#### 4. FASTQC for downloaded file
#### 5. [BWA] Alignment 
#### 6. Generate HOMER's TagDirectory
#### 7. Mapping stats recollection:

The reason this particular pipeline is so simple is because after the HOMER tag directory, the downstream analysis really varies depending on the biological question.

WGB-Seq
-------

_Whole-Genome Bisulfite followed by Sequencing_ is a marvelous and expensive technique that allows the study of methylation analysis at base resolution. It cannot differentiate between 5hmC and 5mC since both almost no not interact with the Sodium Bisulphite treatment the DNA is subjected to.
Therefore, unless OxBS-Seq or TAB-Seq (none of them discussed here... yet) or another similar technique that allows the base-resolution identification of 5hmC (or 5mC) is used, we have to refer to the WGBS signal as __5mC + 5hmC__.


Analysis overview:

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

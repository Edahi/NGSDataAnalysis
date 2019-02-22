# NGSDataAnalysis
I have several scripts to download and process NSG datasets. It would be a good idea to share it with the world.

ATAC-Seq
--------
Assay for Transposase-Accessible Chromatin followed by Sequencing focus on the obtention of genomic regions that are not being wrapped by either nucleosomes or other cellular elements, meaning that they are "Accessible" for interaction with cellular elements.

### 0- (Optional) Download SRR fastq data files
- Download the compressed FASTQ files using [SRA-toolkit's](https://github.com/ncbi/sra-tools) fastq-dump.
### 1- Quality check
- Sequencing lanes are merged before quality assessment
- General sequencing quality, length, barcodes and different oddities are first analyzed through [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
### 2- Step mapping
_1st-step Mapping_
- _Raw_ reads are aligned to reference genome (e.g. mm10) using [Botwie1](http://bowtie-bio.sourceforge.net/manual.shtml)
- Unmapped reads are processed with [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) through the use of [trim_galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Low-quality bases are removed from 5'/3' and as minimum, 1bp is clipped from the 3'-end.

_2nd-step Mapping_
- _Processed_ reads are mapped to the reference genome with more permissive parametrs, such a higher insertion distance and slight more mismatches rate.
- For comparison, _unmapped_ reads are quality-checked with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Both mapping results are merged and sorted using [samtools](http://samtools.sourceforge.net/).
### 3- Contamination assessment
- Unmapped reads are assembled with [velvet and velveth](https://www.ebi.ac.uk/~zerbino/velvet/) into contigs.
- The longest and more covered _Top 3_ contigs are [_Blast_-ed](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to assess any kind of contamination on the sample.
### 4- Cleaning
- Reads aligning to the mitochodrial chromosome are eliminated through the use of both [samtools](http://samtools.sourceforge.net/) and a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script.
- Reads mapping to [ENCODE's blackliested regions](https://sites.google.com/site/anshulkundaje/projects/blacklists) are removed with [bedtools](https://bedtools.readthedocs.io/en/latest/).
- PCR duplicates are removed with [Picard-tools](https://broadinstitute.github.io/picard/)
### 5- Subcategories obtention
- Subnucleosomal reads are obtained with a cutsom [awk](https://www.gnu.org/software/gawk/manual/gawk.html) script (reads closer than 100bp).
- Tn5 insertion sites are obtained through the use of a perl script
### 6- Peak calling
- Accessible peaks from the subnucleosomal reads data are obtained with [MACS2](https://github.com/taoliu/MACS)

Differential Accessible Regions (DARs).
---------------------------------------
There are multiple ways to detect significant differences among a set of conditions. Here I propose the use of Limmma's normalization methods _voom_. Below the detailed description of how I obtain DARs.
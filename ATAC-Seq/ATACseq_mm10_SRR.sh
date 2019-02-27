#!/bin/bash

# ATAC-seq Standarized pipeline mm10: Version -- 03
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2018.10.22
# Versions logs
# v3.
#  Corrected minor filename conveniences (Add HOMER prefix to HOMER peaks).
#  Corrected missing backslash to make final file consistent with my template format.
#  Added MACS2 program and peaks.
#  Renamed folder "06.HOMER_Peaks" by "06.Peaks" given MACS2 peaks introduction.
# v2.
#  Compressed downloaded FASTQ file
#  Reformat the help message
#  Added summary statistics to track mapping results
#  Changed the order between removing duplicates and filtering blacklisted regions
#  Updated the filename to include the current version of the script
#  Added option to remove temporal files, such as intermediate mapping results
#  Renamed TagAlign to Fragment Length Estimate
#  Added code necessary to generate BigWigs
#  Added option to remove intermediate files at the end
#  Add Genome Browser tracks || Experimental with Tn5_9bp
#  Added code to check if SRR was given

# >>> How to use the program:
#    ATACseq_mm10_SRR_v3 [-c SRR code] [-n Name def:ATAC] [-d DownloadDir def:/BioScratch/edahi] [-a AnalysisDir def:/BioScratch/edahi] [-r Run created job def:no] [-q rao-exclusive queue def:no] [-h help]"

usage()
{
    printf "\nusage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_SRR_v3 [-c SRR code]"
    printf "\n\t[-n | --name      Name                       def:ATAC  ]"
    printf "\n\t[-d | --download  DownloadPath               def:/mnt/beegfs ]"
    printf "\n\t[-a | --analysis  AnalysisPath               def:/mnt/BioAdHoc/Groups/RaoLab/temp ]"
    printf "\n\t[-q | --queue     set 'rao-exclusive' queue  def:'default' ]"
    printf "\n\t[-r | --run       Run created job            def:no   ]"
    printf "\n\t[-k | --keep      Keep intermediate results  def:no   ]"
    printf "\n\t[-h | --help      Show this message and exit ]\n\n"
}


# >>> Check if arguments given:
if [ "$1" == "" ]; then
  usage
  exit 1
fi

# >>> Declare Variables
download=/mnt/beegfs
analysis=/mnt/BioAdHoc/Groups/RaoLab/temp
srr=
name=ATAC
run=0
raoqueue=0
keep=0

# >>> Assign arguments to variables
while [ "$1" != "" ]; do
    case $1 in
        -c | --SRR )            shift
                                srr=$1
                                ;;
        -n | --name )           shift
                                name=$1
                                ;;
        -d | --download )       shift
                                download=$1
                                ;;
        -a | --analysis )       shift
                                analysis=$1
                                ;;
        -q | --queue  )         raoqueue=1
                                ;;
        -k | --keep  )          keep=1
                                ;;
        -r | --run    )         run=1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# >>> Check if SRR was given:
if [ "$srr" = "" ]; then
    printf "\n\tSRR CODE is a required argument. Check --help for further assistance\n\n";
    exit
fi

      Jobs=$analysis/Jobs
  softlink=$analysis/01.Data
   mapping=$analysis/02.Mapping
   TagDirs=$analysis/03.TagDirectories
FragLenEst=$analysis/04.FragmentLengthEstimates
       ssp=$analysis/05.SSP
   peakDir=$analysis/06.Peaks
   BigWigs=$analysis/07.BigWigs
    FASTQC=$analysis/FASTQC
  FASTQCun=$analysis/FASTQC_UnMap

mkdir -p $Jobs

cat <<EOF> $Jobs/${name}.ATACseq.mm10.v3.sh
#!/bin/bash -ex
#PBS -N ${name}.ATACseq.mm10.v3
#PBS -l walltime=168:00:00
#PBS -o $Jobs/${name}.ATACseq.mm10.v3.out
#PBS -e $Jobs/${name}.ATACseq.mm10.v3.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
EOF
if [ "$raoqueue" = "1" ]; then
 cat <<EOF>> $Jobs/${name}.ATACseq.mm10.v3.sh
#PBS -q rao-exclusive
EOF
else
 cat <<EOF>> $Jobs/${name}.ATACseq.mm10.v3.sh
#PBS -q default
EOF
fi

cat <<EOF>> $Jobs/${name}.ATACseq.mm10.v3.sh

export PATH=/share/apps/R/3.1.0/bin:/share/apps/python/python-3.4.6/bin:/share/apps/python/python-2.7.13/bin:/share/apps/perl/perl-5.18.1-threaded/bin/:/share/apps/gcc/6.3/bin:/mnt/BioApps/pigz/latest/bin:/share/apps/bin:/usr/local/maui/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/stack/bin:/share/apps/java/latest/bin:/share/apps/bowtie/bowtie2-2.1.0:/share/apps/bowtie/bowtie-1.1.2:/usr/local/cuda/bin:/share/apps/dos2unix/latest/bin:/share/apps/bedtools/bin:/share/apps/HOMER/bin
printf "PATH Used:\n\$PATH\n\n"
unset PYTHONPATH

#Variables:
    bowtie=/Bioinformatics/apps/bowtie/bowtie-1.0.0/bowtie
    fastqc=/Bioinformatics/apps/FastQC/FastQC-0.11.2/fastqc
   MarkDup=/Bioinformatics/apps/picard-tools/picard-tools-1.94/MarkDuplicates.jar
    BamCov=/home/edahi/.local/bin/bamCoverage
       Ann=/home/edahi/download/code/HOMER/bin/annotatePeaks.pl
  getPeaks=/home/edahi/download/code/HOMER/bin/findPeaks
 mMultiWig=/home/edahi/download/code/HOMER/bin/makeMultiWigHub.pl
   makeTag=/home/edahi/download/code/HOMER/bin/makeTagDirectory
 fastqdump=/home/edahi/download/code/sra/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump
trimgalore=/home/edahi/download/code/TrimGalore/0.3.8/trim_galore
   mm10bwa=/home/edahi/download/genome/mm10/bwa/index/mm10_random/mm10.fa
mm10genome=/mnt/BioAdHoc/Groups/RaoLab/Bioinformatics/apps/Mus_musculus/UCSC/mm10/BowtieIndex/genome
genomesize=/mnt/BioAdHoc/Groups/RaoLab/Bioinformatics/apps/Mus_musculus/UCSC/mm10/BowtieIndex/genome.fa.sizes
    Btools=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/bedtools2/bin/bedtools
    makeBG=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/HOMER/bin/makeUCSCfile
    RunSPP=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/phantompeakqualtools/run_spp.R
    BLmm10=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/Mus_musculus/UCSC/mm10/Sequence/Blacklist/mm10.blacklist.bed
      py27=/share/apps/python/python-2.7.13/bin/python
     MACS2=/share/apps/python/python-2.7.6/bin/macs2
   Rscript=/share/apps/R/3.1.0/bin/Rscript
     bg2bw=/share/apps/UCSC/bedGraphToBigWig
      calc=/share/apps/UCSC/calc
       sam=/usr/bin/samtools


# From Script:
  download=$download
       srr=$srr
      name=$name
  softlink=$softlink
   mapping=$mapping
   TagDirs=$TagDirs
FragLenEst=$FragLenEst
       ssp=$ssp
   peakDir=$peakDir
   BigWigs=$BigWigs
    FASTQC=$FASTQC
  FASTQCun=$FASTQCun
      keep=$keep
      Jobs=$Jobs

# Generate Directories Framework:
mkdir -p \$softlink \$mapping/\${name} \$TagDirs \$FragLenEst \$ssp \$peakDir \$FASTQC \$FASTQCun \$BigWigs

# >>>>>> Check if files existed before (THEN DELETOS)
if [ -f \$download/\${srr}_1.fastq.gz ];     then printf "Removing Pre-existing \$download/\${srr}_[12].fastq.gz"     ; rm \$download/\${srr}_*.fastq.gz    ; fi
if [ -L \${softlink}/\${name}_R1.fastq.gz ]; then printf "Removing Pre-existing \${softlink}/\${name}_R[12].fastq.gz" ; rm \${softlink}/\${name}_R*.fastq.gz; fi

# >>>>>> Download Files for ${name} -- $srr
\$fastqdump --dumpbase --skip-technical --clip --outdir \$download --gzip -A \$srr --split-3  --split-files

# >>>>>> Create Softlinks for ${name}
ln -s \$download/${srr}_1.fastq.gz \${softlink}/${name}_R1.fastq.gz
ln -s \$download/${srr}_2.fastq.gz \${softlink}/${name}_R2.fastq.gz

# >>>>>> FASTQC for downloaded files
\$fastqc \${softlink}/\${name}_R1.fastq.gz \${softlink}/\${name}_R2.fastq.gz --outdir=\$FASTQC &

# >>>>>> Mapping (1)
cd       \$mapping/\${name}
printf "Statistics for bowtie mapping of untrimmed reads \n" > \$mapping/\$name/BowtieStats.txt
nohup \$bowtie -p 4 -m 1 --best --strata -X 2000 \\
  -S --fr --chunkmbs 1024 \$mm10genome \\
  -1  <(zcat \${softlink}/\${name}_R1.fastq.gz) \\
  -2  <(zcat \${softlink}/\${name}_R2.fastq.gz) \\
  \$mapping/\$name/\${name}_mm10.sam \\
  --un \$mapping/\$name/\${name}_unmapped.fastq &>> \$mapping/\$name/BowtieStats.txt

# >>>>>> Summary Statistics (0)(1) -- Name & Total reads
echo \${name} > \$mapping/\${name}/SummaryStats00.txt
grep processed \$mapping/\$name/BowtieStats.txt | cut -f2 -d: | tr -d ' '  > \$mapping/\${name}/SummaryStats01.txt

# >>>>>> trim_galore unmapped reads
\$trimgalore --paired --nextera --length 37 \\
  --stringency 3 --three_prime_clip_R1 1 --three_prime_clip_R2 1 \\
  \$mapping/\$name/\${name}_unmapped_1.fastq \$mapping/\$name/\${name}_unmapped_2.fastq -o \$mapping/\$name

# >>>>>> FASTQC for Unmapped reads
\$fastqc \$mapping/\$name/\${name}_unmapped_1.fastq    \$mapping/\$name/\${name}_unmapped_2.fastq    --outdir=\$FASTQCun &

# >>>>>> Remap filtered-unmapped reads
printf "Statistics for bowtie mapping of trim_galore unmapped reads \n" >> \$mapping/\$name/BowtieStats.txt
nohup \$bowtie -p 4 -m 1 --best --strata -X 2000 \\
  -S --fr --chunkmbs 1024 \$mm10genome \\
  -1 \$mapping/\$name/\${name}_unmapped_1_val_1.fq \\
  -2 \$mapping/\$name/\${name}_unmapped_2_val_2.fq \\
  \$mapping/\$name/\${name}_remapTrimUnmapped_mm10.sam &>> \$mapping/\$name/BowtieStats.txt

# >>>>>> FASTQC for Trim_Galore filtered reads
\$fastqc \$mapping/\$name/\${name}_unmapped_1_val_1.fq \$mapping/\$name/\${name}_unmapped_2_val_2.fq --outdir=\$FASTQCun &

# In this case the warning "[bam_header_read] EOF marker is absent. The input is probably truncated."
# Should be ignored. Using samtools view on pipelines gives uncompressed bam, and these do not have the EOF marker
# >>>>>> Generating 1st step bam
\$sam view -bS \$mapping/\$name/\${name}_mm10.sam | \\
  \$sam view -h -F 4 -b - | \\
  \$sam sort - \$mapping/\$name/\${name}_mm10_onlymapped_sorted

# >>>>>> Generating 2nd step bam (remapped trimgalore)
\$sam view -bS \$mapping/\$name/\${name}_remapTrimUnmapped_mm10.sam | \\
  \$sam view -h -F 4 -b - | \\
  \$sam sort - \$mapping/\$name/\${name}_remapTrimUnmapped_mm10_onlymapped_sorted

# >>>>>> Merge both mapping results
\$sam merge \\
  \$mapping/\${name}/\${name}.mapped.sorted.merged.bam \\
  \$mapping/\$name/\${name}_mm10_onlymapped_sorted.bam \\
  \$mapping/\$name/\${name}_remapTrimUnmapped_mm10_onlymapped_sorted.bam

# >>>>>> Summary Statistics (2) -- First Mapped Reads 
# grep Reported \$mapping/\$name/stats1.txt | cut -f2 -d\\   > \$mapping/\${name}/SummaryStats02.txt
\$sam view -c \$mapping/\$name/\${name}_mm10_onlymapped_sorted.bam   > \$mapping/\${name}/SummaryStats02.txt

# >>>>>> Summary Statistics (3) -- Second Mapped Reads 
# grep Reported \$mapping/\$name/stats2.txt | cut -f2 -d\\   > \$mapping/\${name}/SummaryStats03.txt
\$sam view -c \$mapping/\$name/\${name}_remapTrimUnmapped_mm10_onlymapped_sorted.bam   > \$mapping/\${name}/SummaryStats03.txt

# >>>>>> Summary Statistics (4) -- Total Merged Mapped Reads
\$sam view -c \$mapping/\${name}/\${name}.mapped.sorted.merged.bam > \$mapping/\${name}/SummaryStats04.txt

# >>>>>> FragLenEst (Raw):
  \$sam view -F 0x0204 -o - \$mapping/\${name}/\${name}.mapped.sorted.merged.bam | 
    awk 'BEGIN{OFS="\t"}{if (and(\$2,16) > 0) {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","-"} else {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","+"} }' | 
    gzip -c > \$FragLenEst/\${name}.Raw.tagAlign.gz 

# >>>>>> Run ssp [PhantomPeaks] (Raw)
\$Rscript \$RunSPP \\
    -s=-100:5:600 \\
    -c=\$FragLenEst/\${name}.Raw.tagAlign.gz -savp \\
    -out=\$ssp/\${name}.Raw.FragLenEst.tab 

# >>>>>> fragmentLengthEstimate (Raw)
cat \$ssp/\${name}.Raw.FragLenEst.tab |awk '{print \$3}'|cut -d ',' -f 1 > \$FragLenEst/\${name}.Raw.FragLenEst.cat

# >>>>>> Remove chrM
\$sam view -h \$mapping/\${name}/\${name}.mapped.sorted.merged.bam | \\
  perl -lane 'print \$_ if \$F[2] ne "chrM"' | \\
  \$sam view -bS - > \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.bam

# >>>>>> Summary Statistics (5) -- Mapped reads w/o chrM
\$sam view -c \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.bam  > \$mapping/\${name}/SummaryStats05.txt

# >>>>>> Remove Duplicates:
java -jar \$MarkDup \
   INPUT=\$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.bam \\
   OUTPUT=\$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.rmdup.bam \\
   METRICS_FILE=\$mapping/\${name}/\${name}_PicardMetrics.txt \\
   REMOVE_DUPLICATES=true \\
   ASSUME_SORTED=true

# >>>>>> Summary Statistics (6) -- Deduplicated mapped reads w/o chrM
\$sam view -c \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.rmdup.bam > \$mapping/\${name}/SummaryStats06.txt

# >>>>>> Whitelist:
\$Btools intersect -a \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.rmdup.bam -b \$BLmm10 -v > \$mapping/\${name}/\${name}.whitelist.bam

# >>>>>> Summary Statistics (7) -- Raw clean mapping results
\$sam view -c \$mapping/\${name}/\${name}.whitelist.bam  > \$mapping/\${name}/SummaryStats07.txt

# >>>>>> FragLenEst (Filtered):
  \$sam view -F 0x0204 -o - \$mapping/\${name}/\${name}.whitelist.bam | 
    awk 'BEGIN{OFS="\t"}{if (and(\$2,16) > 0) {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","-"} else {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","+"} }' | 
    gzip -c > \$FragLenEst/\${name}.Filtered.tagAlign.gz 

# >>>>>> Run ssp [PhantomPeaks] (Filtered)
\$Rscript \$RunSPP \\
    -s=-100:5:600 \\
    -c=\$FragLenEst/\${name}.Filtered.tagAlign.gz -savp \\
    -out=\$ssp/\${name}.Filtered.FragLenEst.tab 

# >>>>>> fragmentLengthEstimate (Filtered)
cat \$ssp/\${name}.Filtered.FragLenEst.tab |awk '{print \$3}'|cut -d ',' -f 1 > \$FragLenEst/\${name}.Filtered.FragLenEst.cat

# >>>>>> Extract sub-nucleosomal fragments'
\$sam view -H \$mapping/\${name}/\${name}.whitelist.bam > \$mapping/\${name}/\${name}.subNuc.sam

\$sam view \$mapping/\${name}/\${name}.whitelist.bam | \\
  awk '{if(sqrt(\$9*\$9)<100){print \$0}}' >> \$mapping/\${name}/\${name}.subNuc.sam

\$sam view -S -b \$mapping/\${name}/\${name}.subNuc.sam > \$mapping/\${name}/\${name}.subNuc.bam

# >>>>>> Summary Statistics (8) -- Subnucleosomal Fragments
\$sam view -c \$mapping/\${name}/\${name}.subNuc.bam  > \$mapping/\${name}/SummaryStats08.txt

# >>>>>> Obtain Tn5 footprint (1)BamToBed (2)perlBed (3)BedToBam
\$Btools bamtobed -i \$mapping/\${name}/\${name}.whitelist.bam > \$mapping/\${name}/\${name}.whitelist.bed

/home/edahi/download/code/ATACseq/Tn5_bed9bp_full.pl \$mapping/\${name}/\${name}.whitelist.bed \$mapping/\${name}/\${name}_Tn5footprint_unsort.bed

\$Btools bedtobam -i \$mapping/\${name}/\${name}_Tn5footprint_unsort.bed -g \$genomesize | \\
  \$sam sort - \$mapping/\$name/\${name}.Tn5footprint

# >>>>>> Summary Statistics (9) -- Tn5 (9bp) Insertion sites
\$sam view -c \$mapping/\${name}/\${name}.Tn5footprint.bam  > \$mapping/\${name}/SummaryStats09.txt

# >>>>>> Gather all summaries:
echo 'SampleName,RawReads,Mapped_1st,Mapped_2nd,TotalMapped,chrM_Filter,Duplicated_Filter,Whitelist_Filter,Subnucleosomal,Tn5_9bp_InsertionSizes' >  \$mapping/\${name}/Colnames.MappingStats.csv
paste -d, \$mapping/\${name}/SummaryStats0[0-9].txt  >  \$mapping/\${name}/MappingStats.csv

# >>>>>>Indexes of all files (1)Clean Mapping Results (2)SubNucleosomal (3)Tn5Footprint
\$sam index \$mapping/\${name}/\${name}.whitelist.bam
\$sam index \$mapping/\${name}/\${name}.subNuc.bam
\$sam index \$mapping/\${name}/\${name}.Tn5footprint.bam

# >>>>>> Calculate fragment length distributions
\$py27 /home/edahi/download/code/ATACseq/Fragment_length_density_plot.py \$mapping/\${name}/\${name}.whitelist.bam \$name \$FragLenEst/\${name}_fragmentLengths

# >>>>>> Generate Tag Directories
mkdir -p  \$TagDirs/\${name}_whitelist \$TagDirs/\${name}_subNuc \$TagDirs/\${name}_Tn5
\$makeTag \$TagDirs/\${name}_whitelist -keepAll -illuminaPE \$mapping/\${name}/\${name}.whitelist.bam    > \$TagDirs/\${name}_whitelist/maketagdir.txt
\$makeTag \$TagDirs/\${name}_subNuc    -keepAll -illuminaPE \$mapping/\${name}/\${name}.subNuc.bam       > \$TagDirs/\${name}_subNuc/maketagdir.txt
\$makeTag \$TagDirs/\${name}_Tn5       -keepAll -illuminaPE \$mapping/\${name}/\${name}.Tn5footprint.bam > \$TagDirs/\${name}_Tn5/maketagdir.txt

# >>>>>> Call Peaks
\$getPeaks \$TagDirs/\${name}_whitelist -style dnase -region -nfr -o \$peakDir/HOMER.\${name}.whitelist.peaks.txt
\$getPeaks \$TagDirs/\${name}_subNuc    -style dnase -region -nfr -o \$peakDir/HOMER.\${name}.subNuc.peaks.txt
\$getPeaks \$TagDirs/\${name}_Tn5       -style dnase -region -nfr -o \$peakDir/HOMER.\${name}.Tn5.peaks.txt

\$MACS2 callpeak -t \$mapping/\${name}/\${name}.whitelist.bam    -f BAMPE  -n \$peakDir/MACS2.\${name}.whitelist -g mm -q 0.05  --keep-dup all --nomodel
\$MACS2 callpeak -t \$mapping/\${name}/\${name}.subNuc.bam       -f BAMPE  -n \$peakDir/MACS2.\${name}.subNuc    -g mm -q 0.05  --keep-dup all --nomodel
\$MACS2 callpeak -t \$mapping/\${name}/\${name}.Tn5footprint.bam -f BAM    -n \$peakDir/MACS2.\${name}.Tn5       -g mm -q 0.05  --keep-dup all --nomodel

# >>>>>> Generate BigWig files 
\$makeBG   \$TagDirs/\${name}_whitelist -o auto -fsize 1e20 
\$makeBG   \$TagDirs/\${name}_subNuc    -o auto -fsize 1e20 
\$makeBG   \$TagDirs/\${name}_Tn5       -o auto -fsize 1e20 

zcat \$TagDirs/\${name}_whitelist/\${name}_whitelist.ucsc.bedGraph.gz | tail -n+2 | sort -k1,1 -k2,2n > \$TagDirs/\${name}_whitelist/\${name}_whitelist.ucsc.bedGraph
zcat \$TagDirs/\${name}_subNuc/\${name}_subNuc.ucsc.bedGraph.gz       | tail -n+2 | sort -k1,1 -k2,2n > \$TagDirs/\${name}_subNuc/\${name}_subNuc.ucsc.bedGraph
zcat \$TagDirs/\${name}_Tn5/\${name}_Tn5.ucsc.bedGraph.gz             | tail -n+2 | sort -k1,1 -k2,2n > \$TagDirs/\${name}_Tn5/\${name}_Tn5.ucsc.bedGraph

\$bg2bw \$TagDirs/\${name}_whitelist/\${name}_whitelist.ucsc.bedGraph \$genomesize \$BigWigs/\${name}.whitelist.bw
\$bg2bw \$TagDirs/\${name}_subNuc/\${name}_subNuc.ucsc.bedGraph       \$genomesize \$BigWigs/\${name}.subNuc.bw
\$bg2bw \$TagDirs/\${name}_Tn5/\${name}_Tn5.ucsc.bedGraph             \$genomesize \$BigWigs/\${name}.Tn5.bw

# \$BamCov --bam \$mapping/\${name}/\${name}.Tn5footprint.bam --numberOfProcessors 4 --binSize 1  --normalizeUsing RPKM --smoothLength 1  --ignoreDuplicates -o \$BigWigs/\${name}.Tn5_9bp.bw  --maxFragmentLength 10 &
# \$BamCov --bam \$mapping/\${name}/\${name}.subNuc.bam       --numberOfProcessors 4 --binSize 10 --normalizeUsing RPKM --smoothLength 30 --ignoreDuplicates -o \$BigWigs/\${name}.subNuc.bw &
# \$BamCov --bam \$mapping/\${name}/\${name}.whitelist.bam    --numberOfProcessors 4 --binSize 10 --normalizeUsing RPKM --smoothLength 30 --ignoreDuplicates -o \$BigWigs/\${name}.whitelist.bw 

# >>> Add BigWigs to tracks file:
echo track type=bigWig name=\${name}.whitelist description=\${name}.whitelist visibility=2 autoScale=off maxHeightPixels=40 viewLimits=0:200 color=10,10,10 graphType=bar bigDataUrl=http://informaticsdata.liai.org/NGS_analyses/ad_hoc/\${BigWigs/\\/mnt\/BioAdHoc\\//}/\${name}.whitelist.bw >> \${BigWigs}/ATAC_Tracks_Whitelist.txt
echo track type=bigWig name=\${name}.subNuc    description=\${name}.subNuc    visibility=2 autoScale=off maxHeightPixels=40 viewLimits=0:200 color=0,0,255  graphType=bar bigDataUrl=http://informaticsdata.liai.org/NGS_analyses/ad_hoc/\${BigWigs/\\/mnt\/BioAdHoc\\//}/\${name}.subNuc.bw    >> \${BigWigs}/ATAC_Tracks_subNuc.txt
echo track type=bigWig name=\${name}.Tn5       description=\${name}.Tn5       visibility=2 autoScale=off maxHeightPixels=40 viewLimits=0:200 color=0,255,0  graphType=bar bigDataUrl=http://informaticsdata.liai.org/NGS_analyses/ad_hoc/\${BigWigs/\\/mnt\/BioAdHoc\\//}/\${name}.Tn5_9bp.bw   >> \${BigWigs}/ATAC_Tracks_Tn5_9bp.txt

# >>>>>> Remove Intermediate Files (if selected):
if [ "\$keep" = "0" ]; then
 rm \$mapping/\$name/\${name}*.sam
 rm \$mapping/\$name/\${name}_unmapped*q
 rm \$mapping/\${name}/\${name}*mapped*.bam
 rm \$FragLenEst/\${name}.*.tagAlign.gz 
 rm \$mapping/\${name}/\${name}*.bed
 rm \$mapping/\${name}/SummaryStats0[0-9].txt
 rm \$mapping/\${name}/gmon.out
fi

EOF
  
# >>> If selected, Run job
if [ "$run" = "1" ]; then
 echo "Running Job"
 qsub $Jobs/${name}.ATACseq.mm10.v3.sh
fi

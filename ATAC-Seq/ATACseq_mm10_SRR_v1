#!/bin/bash

# ATAC-seq Standarized pipeline mm10: Version -- 01
# >>> How to use the program:
#    ATACseq_mm10_SRR_v1 [-c SRR code] [-n Name def:ATAC] [-d DownloadDir def:/BioScratch/edahi] [-a AnalysisDir def:/BioScratch/edahi] [-r Run created job def:no] [-h help]"

usage()
{
    echo "usage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_SRR_v1 [-c SRR code] [-n Name def:ATAC] [-d DownloadDir def:/BioScratch/edahi]  [-a AnalysisDir def:/BioScratch/edahi]  [-r Run created job def:no] [-h help]"
}

# >>> Check if arguments given:
if [ "$1" == "" ]; then
  usage
  exit 1
fi

# >>> Declare Variables
download=/BioScratch/edahi
analysis=/BioScratch/edahi
srr=
name=ATAC
run=0

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

    Jobs=$analysis/Jobs
softlink=$analysis/01.Data
 mapping=$analysis/02.Mapping
 TagDirs=$analysis/03.TagDirectories
TagAlign=$analysis/04.TagAlign
     ssp=$analysis/05.SSP
 peakDir=$analysis/06.HOMER_Peaks
  FASTQC=$analysis/FASTQC

mkdir -p $Jobs

cat <<EOF> $Jobs/${name}.ATACseq.mm10.sh
#!/bin/bash -ex
#PBS -N ${name}.ATACseq.mm10
#PBS -l walltime=168:00:00
#PBS -o $Jobs/${name}.ATACseq.mm10.out
#PBS -e $Jobs/${name}.ATACseq.mm10.out
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -q rao-exclusive
#PBS -m ae

export PATH=/share/apps/R/3.1.0/bin:/share/apps/python/python-3.4.6/bin:/share/apps/python/python-2.7.13/bin:/share/apps/perl/perl-5.18.1-threaded/bin/:/share/apps/gcc/6.3/bin:/mnt/BioApps/pigz/latest/bin:/share/apps/bin:/usr/local/maui/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/stack/bin:/share/apps/java/latest/bin:/share/apps/bowtie/bowtie2-2.1.0:/share/apps/bowtie/bowtie-1.1.2:/usr/local/cuda/bin:/share/apps/dos2unix/latest/bin:/share/apps/bedtools/bin:/share/apps/HOMER/bin
echo \$PATH
unset PYTHONPATH

#Variables:
    bowtie=/Bioinformatics/apps/bowtie/bowtie-1.0.0/bowtie
    fastqc=/Bioinformatics/apps/FastQC/FastQC-0.11.2/fastqc
   MarkDup=/Bioinformatics/apps/picard-tools/picard-tools-1.94/MarkDuplicates.jar
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
    BLmm10=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/Mus_musculus/UCSC/mm10/Sequence/Blacklist/mm10.blacklist.bed
       sam=/usr/bin/samtools

# From Script:
  download=$download
       srr=$srr
      name=$name
  softlink=$softlink
   mapping=$mapping
   TagDirs=$TagDirs
  TagAlign=$TagAlign
       ssp=$ssp
   peakDir=$peakDir
    FASTQC=$FASTQC
      Jobs=$Jobs

# Generate Directories Framework:
mkdir -p \$softlink \$mapping/\${name} \$TagDirs \$TagAlign \$ssp \$peakDir \$FASTQC

# >>>>>> Check if files existed before (THEN DELETOS)
if [ -f \$download/\${srr}_1.fastq ];    then rm \$download/\${srr}_*.fastq   ; fi
if [ -f \${softlink}/\${name}_R1.fastq ]; then rm \${softlink}/\${name}_R*.fastq; fi

# >>>>>> Download Files for ${name} -- $srr
\$fastqdump --dumpbase --split-files --skip-technical --clip --outdir \$download --split-3 -A \$srr

# >>>>>> Create Softlinks for ${name}
ln -s \$download/${srr}_1.fastq \${softlink}/${name}_R1.fastq
ln -s \$download/${srr}_2.fastq \${softlink}/${name}_R2.fastq

# >>>>>> FASTQC for downloaded files
\$fastqc \${softlink}/\${name}_R1.fastq \${softlink}/\${name}_R2.fastq --outdir=\$FASTQC

# >>>>>> Mapping (1)
cd       \$mapping/\${name}
printf "Statistics for bowtie mapping of untrimmed reads \n" > \$mapping/\$name/stats.txt
nohup \$bowtie -p 4 -m 1 --best --strata -X 2000 \\
  -S --fr --chunkmbs 1024 \$mm10genome \\
  -1  \${softlink}/\${name}_R1.fastq \\
  -2  \${softlink}/\${name}_R2.fastq \\
  \$mapping/\$name/\${name}_mm10.sam \\
  --un \$mapping/\$name/\${name}_unmapped.fastq &>> \$mapping/\$name/stats.txt

# >>>>>> trim_galore unmapped reads
\$trimgalore --paired --nextera --length 37 \\
  --stringency 3 --three_prime_clip_R1 1 --three_prime_clip_R2 1 \\
  \$mapping/\$name/\${name}_unmapped_1.fastq \$mapping/\$name/\${name}_unmapped_2.fastq -o \$mapping/\$name

# >>>>>> Remap filtered-unmapped reads
printf "Statistics for bowtie mapping of trim_galore unmapped reads \n" >> \$mapping/\$name/stats.txt
nohup \$bowtie -p 4 -m 1 --best --strata -X 2000 \\
  -S --fr --chunkmbs 1024 \$mm10genome \\
  -1 \$mapping/\$name/\${name}_unmapped_1_val_1.fq \\
  -2 \$mapping/\$name/\${name}_unmapped_2_val_2.fq \\
  \$mapping/\$name/\${name}_remapTrimUnmapped_mm10.sam &>> \$mapping/\$name/stats.txt

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

# >>>>>> Remove chrM
\$sam view -h \$mapping/\${name}/\${name}.mapped.sorted.merged.bam | \\
  perl -lane 'print \$_ if \$F[2] ne "chrM"' | \\
  \$sam view -bS - > \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.Blacklist.bam

# >>>>>> Whitelist:
\$Btools intersect -a \$mapping/\${name}/\${name}.mapped.sorted.merged.nochrM.Blacklist.bam -b \$BLmm10 -v > \$mapping/\${name}/\${name}.whitelist.bam

# >>>>>> Remove Duplicates:
java -jar \$MarkDup \
   INPUT=\$mapping/\${name}/\${name}.whitelist.bam \
   OUTPUT=\$mapping/\${name}/\${name}.whitelist_rmdup.bam \
   METRICS_FILE=\$mapping/\${name}/\${name}_PicardMetrics.txt \
   REMOVE_DUPLICATES=true \
   ASSUME_SORTED=true

# >>>>>> TagAlign:
  \$sam view -F 0x0204 -o - \$mapping/\${name}/\${name}.whitelist_rmdup.bam | 
    awk 'BEGIN{OFS="\t"}{if (and(\$2,16) > 0) {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","-"} else {print \$3,(\$4-1),(\$4-1+length(\$10)),"N","1000","+"} }' | 
    gzip -c > \$TagAlign/\${name}.mapped.tagAlign.gz 

# >>>>>> Run ssp
/share/apps/R/3.1.0/bin/Rscript /home/edahi/download/code/phantompeakqualtools/run_spp.R \
    -s=-100:5:600 \
    -c=\$TagAlign/\${name}.mapped.tagAlign.gz -savp \
    -out=\$ssp/\${name}.mapped.tagAlign.tab 

# >>>>>> UndetCat
cat \$TagAlign/\${name}.mapped.tagAlign.tab |awk '{print \$3}'|cut -d ',' -f 1 > \$TagAlign/\${name}.mapped.tagAlign.cat

# >>>>>> Extract sub-nucleosomal fragments'
\$sam view -H \$mapping/\${name}/\${name}.whitelist_rmdup.bam > \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.sam

\$sam view \$mapping/\${name}/\${name}.whitelist_rmdup.bam | \\
  awk '{if(sqrt(\$9*\$9)<100){print \$0}}' >> \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.sam

\$sam view -S -b \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.sam > \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam

# >>>>>> Obtain Tn5 footprint (1)BamToBed (2)perlBed (3)BedToBam
\$Btools bamtobed -i \$mapping/\${name}/\${name}.whitelist_rmdup.bam > \$mapping/\${name}/\${name}.whitelist_rmdup.bed

# >>>>>> Isolate Tn5footprint'
/home/edahi/download/code/ATACseq/Tn5_bed9bp_full.pl \$mapping/\${name}/\${name}.whitelist_rmdup.bed \$mapping/\${name}/\${name}_mm10_merge_onlymapped_sorted_rmdup_nochrM_Tn5footprint.bed

# >>>>>> Generate bam from Tn5footprint BED'
\$Btools bedtobam -i \$mapping/\${name}/\${name}_mm10_merge_onlymapped_sorted_rmdup_nochrM_Tn5footprint.bed -g \$genomesize | \\
  \$sam sort - \$mapping/\$name/\${name}_mm10_merge_onlymapped_sorted_rmdup_nochrM_Tn5footprint_sorted

# Indexes of all files (1)Clean Mapping Results (2)SubNucleosomal (3)Tn5Footprint
# >>>>>> Index from Clean Mapping Results
\$sam index \$mapping/\${name}/\${name}.whitelist_rmdup.bam
# >>>>>> Index from subNuc fragments
\$sam index \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam
# >>>>>> Index from sorted 9bp footprint bam
\$sam index \$mapping/\$name/\${name}_mm10_merge_onlymapped_sorted_rmdup_nochrM_Tn5footprint_sorted.bam

# >>>>>> Calculate fragment length distributions
python /home/edahi/download/code/ATACseq/Fragment_length_density_plot.py \$mapping/\${name}/\${name}.whitelist_rmdup.bam \$name \$mapping/\${name}/\${name}_fragmentLengths

# >>>>>> subNuc HOMER tagDirectory
mkdir \$TagDirs/\${name}_subNuc
cd    \$TagDirs/\${name}_subNuc
\$makeTag \$TagDirs/\${name}_subNuc -keepAll -illuminaPE \$mapping/\${name}/\${name}_subNuc_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam > \$TagDirs/\${name}_subNuc/maketagdir.txt

# >>>>>> Tn5 HOMER tagDirectory
mkdir \$TagDirs/\${name}_Tn5tagDir
cd    \$TagDirs/\${name}_Tn5tagDir
\$makeTag \$TagDirs/\${name}_Tn5tagDir -keepAll -illuminaPE \$mapping/\$name/\${name}_mm10_merge_onlymapped_sorted_rmdup_nochrM_Tn5footprint_sorted.bam > \$TagDirs/\${name}_Tn5tagDir/maketagdir.txt

# >>>>>> Raw HOMER tagDirectory
mkdir \$TagDirs/\${name}_Raw
cd    \$TagDirs/\${name}_Raw
\$makeTag \$TagDirs/\${name}_Raw -keepAll -illuminaPE \$mapping/\${name}/\${name}.whitelist_rmdup.bam > \$TagDirs/\${name}_Raw/maketagdir.txt

# >>>>>> Peaks from subNuc tags
\$getPeaks \$TagDirs/\${name}_subNuc  -style dnase -region -nfr -o \$peakDir/\${name}_subNuc.peaks.txt

# >>>>>> Peaks from Tn5 tags
\$getPeaks \$TagDirs/\${name}_Tn5tagDir -style dnase -region -nfr -o \$peakDir/\${name}_Tn5tagDir.peaks.txt

# >>>>>> Peaks from Raw tags
\$getPeaks \$TagDirs/\${name}_Raw  -style dnase -region -nfr -o \$peakDir/\${name}_Raw.peaks.txt

EOF

  
# >>> Rscript to calculate Stats:
if [ "$run" = "1" ]; then
 echo "Running Job"
 qsub $Jobs/${name}.ATACseq.mm10.sh
fi

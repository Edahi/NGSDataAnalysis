#!/bin/bash

# ATAC-seq BaseSpace Download Standarized pipeline hg38: Version -- 01
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2018.09.24
# Description: download fastq from Illumina's BaseSpace and analize data.
# Versions logs
# v1.
#  Download data, create scripts for each downloaded sample and run it
#  Current version linked to:
#  1) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_hg38_Existing_v1
#
# >>> How to use the program:
#    ATACseq_hg38_BSdownload_v1 [-b BaseSpace_project_ID || -t Run TroubleShooting] [-d DownloadDir def:/BioScratch/edahi] [-a AnalysisDir def:/BioScratch/edahi] [-r Run created job def:no] [-q rao-exclusive queue def:no] [-h help]"
# >>> Example
# $E/00.Scripts/Bash/ATACseq_hg38_BSdownload_v1.sh -b 94124062 -d /mnt/BioScratch/edahi/Giuliana/ATACseq/03.hg38_NaiveCells_EmptyCANFAT1_and_CARITNFAT1/ -a $E/Projects/Giuliana/ATAC/03.hg38_NaiveCells_EmptyCANFAT1_and_CARITNFAT1/ -q -r
# $E/00.Scripts/Bash/ATACseq_hg38_BSdownload_v1.sh -b 92010932 -d /mnt/BioScratch/edahi/Hiroshi/ATACseq/01.8_22_18_HiroshiATAC -a $E/Projects/Hiroshi/ATAC/01.8_22_18_HiroshiATAC -q -r
# $E/00.Scripts/Bash/ATACseq_hg38_BSdownload_v1.sh -b 43503462 -d /mnt/BioScratch/edahi/Jerry/ATACseq/01.7_14_17_Jerry_ATAC -a $E/Projects/Jerry/ATAC/01.7_14_17_Jerry_ATAC -q -r
# $E/00.Scripts/Bash/ATACseq_hg38_BSdownload_v1.sh -t          -d /mnt/BioScratch/edahi/Giuliana/ATACseq/01.NaiveCells_EmptyCANFAT1_and_CARITNFAT1/ -a $E/Projects/Giuliana/ATAC/01.NaiveCells_EmptyCANFAT1_and_CARITNFAT1/ -q -r -2 -3

usage()
{
    printf "\nusage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_BSdownload_v1 [-b|--BS BaseSpace_project_ID ]"
    printf "\n\t[ -d | --download  DownloadPath                     def:/mnt/beegfs ]"
    printf "\n\t[ -a | --analysis  AnalysisPath                     def:/mnt/BioAdHoc/Groups/RaoLab/temp ]"
    printf "\n\t[ -q | --queue     set 'rao-exclusive' queue        def:'default' ]"
    printf "\n\t[ -r | --run       Run created job                  def:no   ]"
    printf "\n\t[ -k | --keep      Keep intermediate results        def:no   ]"
    printf "\n\t[ -t | --ts        Run TroubleShooting for Library  def:no   ]"
    printf "\n\t[ -1 | --ts1       Set TroubleShooting #1           def:no   ]"
    printf "\n\t[ -2 | --ts2       Set TroubleShooting #2           def:no   ]"
    printf "\n\t[ -3 | --ts3       Set TroubleShooting #3           def:no   ]"
    printf "\n\t[ -h | --help      Show this message and exit                ]\n\n"
}


# >>> Check if arguments given:
if [ "$1" == "" ]; then
  usage
  exit 1
fi

# >>> Declare Variables
temp=/mnt/beegfs/edahi
download=/mnt/beegfs
analysis=/mnt/BioAdHoc/Groups/RaoLab/temp
pid=
name=ATAC
raoqueue=0
run=0
keep=0
ts=0
ts1=
ts2=
ts3=

# >>> Assign arguments to variables
while [ "$1" != "" ]; do
    case $1 in
        -b | --BS )            shift
                                pid=$1
                                ;;
        -d | --download )       shift
                                download=$1
                                ;;
        -a | --analysis )       shift
                                analysis=$1
                                ;;
        -q | --queue  )         raoqueue=1
                                ;;
        -r | --run    )         run=1
                                ;;
        -k | --keep  )          keep=1
                                ;;
        -t | --ts  )            ts=1
                                ;;
        -1 | --ts1  )           ts1="-1"
                                ;;
        -2 | --ts2  )           ts2="-2"
                                ;;
        -3 | --ts3  )           ts3="-3"
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
mkdir -p $Jobs

# >>> Check if TroubleShooting:
if [ "$ts" = "1" ]; then
cat <<EOF> $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.sh
#!/bin/bash -ex
#PBS -N TroubleShooting.ATACseq.BSdownload.hg38.v1
#PBS -l walltime=168:00:00
#PBS -o $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.out
#PBS -e $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
EOF
if [ "$raoqueue" = "1" ]; then
 cat <<EOF>> $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.sh
#PBS -q rao-exclusive
EOF
else
 cat <<EOF>> $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.sh
#PBS -q default
EOF
fi
cat <<EOF>> $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.sh
# From Script:
      pid=$pid
 analysis=$analysis
 download=$download
 raoqueue=$raoqueue
     keep=$keep
     Jobs=$Jobs
     temp=$temp
# Obtain dataset's names:
nameS=(\$(basename -a \$(ls \$download/* | grep L001_R1)))
nameS=("\${nameS[@]/_L001_R1_001.fastq.gz/}")
# Iterate throught the names of the downloaded sets AND
# Build the TroubleShooting script for each file:
for i in \${!nameS[@]}; do 
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g'); # echo \$name -- \${nameS[\$i]}
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_hg38_Existing_TroubleShooting_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis -q -r $ts1 $ts2 $ts3
done
EOF

 if [ "$run" = "1" ]; then
  echo "Running TroubleShooting For Whole Library"
  qsub $Jobs/TroubleShooting.ATACseq.BSdownload.hg38.v1.sh
 fi
 exit
fi

# >>> Check if PID was given:
if [ "$pid" = "" ]; then
    printf "\n\tPID CODE is a required argument if not TroubleShooting. \nCheck --help for further assistance\n\n";
    exit
fi
# >>> Download BaseSpace data:
cat <<EOF> $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.sh
#!/bin/bash -ex
#PBS -N BaseScript.ATACseq.BSdownload.hg38.v1
#PBS -l walltime=168:00:00
#PBS -o $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.out
#PBS -e $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
EOF
if [ "$raoqueue" = "1" ]; then
 cat <<EOF>> $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.sh
#PBS -q rao-exclusive
EOF
else
 cat <<EOF>> $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.sh
#PBS -q default
EOF
fi
cat <<EOF>> $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.sh
# From Script:
      pid=$pid
 analysis=$analysis
 download=$download
 raoqueue=$raoqueue
     keep=$keep
     Jobs=$Jobs
     temp=$temp

mkdir -p \$download
# Download data from BaseSpace Using MY CREDENTIALS (Project must be shared with me in order to work).
perl /share/apps/BSDownload/BSDownload-Generic.pl -c 739a72400375493aac50a8aeb910c3eb -s 7e4fd758cef94f31845a51116f007373 -t 7acae786b29b4a00ae883f6b80edaef5  -o \$download -p \$pid

# Obtain dataset's names:
nameS=(\$(basename -a \$(ls \$download/* | grep L001_R1)))
nameS=("\${nameS[@]/_L001_R1_001.fastq.gz/}")
# name=("\${nameS[@]/_S?([0-9])?([0-9])/}")

# Iterate throught the names of the downloaded sets AND
# Build the analysis script for each file:
for i in \${!nameS[@]}; do 
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g');
 # echo \$name -- \${nameS[\$i]}
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_hg38_Existing_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis -q -r
done
EOF

# >>> If selected, Run job
if [ "$run" = "1" ]; then
 echo "Running Job"
 qsub $Jobs/BaseScript.ATACseq.BSdownload.hg38.v1.sh
fi

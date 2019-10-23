#!/bin/bash

# Whole Genome Sequencing BaseSpace Download Standarized pipeline: Version -- 01
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2019.04.21
# Description: download fastq from Illumina's BaseSpace and analize data.
# Versions logs
# v1.
#  Download data
#  Select genome build
#  Create scripts for each downloaded sample
#
#  Current version linked to:
#  1) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/WGSeq_mm10_PE_Existing_v1
#  2) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPCseq_mm10_PE_Existing_v1
#  3) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPCseq_hg19_SE_Existing_v1
#  4) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPCseq_hg38_PE_Existing_v1
#  5) 
#
# >>> How to use the program:
#    WGSeq_BSdownload_v1.sh [-b BaseSpace_project_ID || -t Run TroubleShooting] [-s Should someone do SE sequencing] [-d DownloadDir def:/BioScratch/edahi] [-a AnalysisDir def:/BioScratch/edahi] [-r Run created job def:no] [-q rao-exclusive queue def:no] [-h help]"
# >>> Example
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/WGSeq_BSdownload_v1.sh -b 126494692 -g mm10 -l PE -d /mnt/BioScratch/edahi/4_19_19_Edahi_Gonzalez -a /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/05.DemethylationDynamics/04.Foxp3_PacBio/03.SequencingTrials/4_19_19_Edahi_Gonzalez -q -r
#

usage()
{
    printf "\nusage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/WGSeq_BSdownload_v1.sh"
    printf "\n\t[ -b | --BS        BaseSpace_project_ID | 'PROXI'                  def:NONE -- Required Argument ]"
    printf "\n\t[ -g | --Genome    Genome Build  mm9|mm10|hg19|hg38                def:NONE -- Required Argument ]"
    printf "\n\t[ -l | --Library   Library Prep  SE|PE  'SingleEnd|PairedEnd'      def:NONE -- Required Argument ]"
    printf "\n\t[ -d | --download  DownloadPath                                    def:/mnt/beegfs ]"
    printf "\n\t[ -a | --analysis  AnalysisPath                                    def:/mnt/BioScratch/edahi/WGSeq ]"
    printf "\n\t[ -q | --queue     set 'rao-exclusive' queue                       def:'default' ]"
    printf "\n\t[ -r | --run       Run created job                                 def:no   ]"
    printf "\n\t[ -k | --keep      Keep intermediate results                       def:no   ]"
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
analysis=/mnt/BioScratch/edahi/WGseq
     pid=
  Genome=
 Library=
    name=WGseq
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
        -b | --BS       )       shift
                                pid=$1
                                ;;
        -g | --Genome   )       shift
                                Genome=$1
                                ;;
        -l | --Library  )       shift
                                Library=$1
                                ;;
        -d | --download )       shift
                                download=$1
                                ;;
        -a | --analysis )       shift
                                analysis=$1
                                ;;
        -q | --queue    )       raoqueue=1
                                ;;
        -r | --run      )       run=1
                                ;;
        -k | --keep     )       keep=1
                                ;;
        -t | --ts       )       ts=1
                                ;;
        -h | --help     )       usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# >>>>>> Other than Troubleshooting:
# >>> Check if PID was given and if Genome build and Library type was correctly set:
Over(){ printf "\n\tBasespace PID is a required argument. \nCheck --help for further assistance\n\n"; exit ;}
[[ $pid == "" ]] && Over   ||  echo '  Basespace project ID: '$pid

list=(mm9 mm10 hg19 hg38)
Over(){ printf "\n\tProvide a valid genome build \n\tCheck --help for further assistance\n\n"; exit ;}
[[ ${list[@]} =~ (^|[[:space:]])$Genome($|[[:space:]]) ]] && echo '          Genome build: '$Genome   ||  Over

list=(SE PE)
Over(){ printf "\n\tProvide a valid Library type \n\tCheck --help for further assistance\n\n"; exit ;}
[[ ${list[@]} =~ (^|[[:space:]])$Library($|[[:space:]]) ]] && echo '          Library type: '$Library ||  Over

# >>> Obtain parameters for script:
Parameters=

if [ "$raoqueue" = "1" ]; then
 Parameters=" "${Parameters}" -q"
fi

if [ "$run" = "1" ]; then
 Parameters=${Parameters}" -r"
fi

if [ "$keep" = "1" ]; then
 Parameters=${Parameters}" -k"
fi

Jobs=$analysis/Jobs
mkdir -p $Jobs

[[ $raoqueue == "1" ]]  && queue=rao-exclusive ||  queue=default

# >>> Download BaseSpace data:
cat <<EOF> $Jobs/BaseScript.WGSeq.BSdownload.v1.sh
#!/bin/bash -ex
#PBS -N BaseScript.WGSeq.BSdownload.v1
#PBS -l walltime=168:00:00
#PBS -o $Jobs/BaseScript.WGSeq.BSdownload.v1.out
#PBS -e $Jobs/BaseScript.WGSeq.BSdownload.v1.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
#PBS -q $queue

# From Script:
      pid=$pid
 analysis=$analysis
 download=$download
 raoqueue=$raoqueue
     keep=$keep
     Jobs=$Jobs
     temp=$temp

# >>>>>> Download data from BaseSpace using MY CREDENTIALS (Project must be shared with me in order to work).
if [ "\$pid" = "PROXI" ]; then
  echo Data downloaded before, skipping download...
else
  echo Downloading \$pid --to-- \$download
  perl /share/apps/BSDownload/BSDownload-Generic.pl -c 739a72400375493aac50a8aeb910c3eb -s 7e4fd758cef94f31845a51116f007373 -t 7acae786b29b4a00ae883f6b80edaef5  -o \$download -p \$pid
fi

# Obtain dataset's names:
nameS=(\$(basename -a \$(ls $download/* | grep L00._R1 | rev | cut -c22- | rev | uniq)))

# Iterate throught the names of the downloaded sets AND Collect qsub codes:
CollectQsubs=

# Build the analysis script for each file:
for i in \${!nameS[@]}; do 
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g');
 echo \$name -- \${nameS[\$i]}
 CollectQsubs=\${CollectQsubs}:\$(/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/WGSeq_${Genome}_${Library}_Existing_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis $Parameters)
done

CollectQsubs=\${CollectQsubs#:}

# >>> After all the jobs are completed:
# >>> Run Collect Mapping Stats
cat <<EOT> $Jobs/SummaryStats.Mapping.v1.sh
#!/bin/bash -x
#PBS -N WGSeq-PostProcessing-v1
#PBS -l walltime=8:00:00
#PBS -o $Jobs/SummaryStats.Mapping.v1.out
#PBS -e $Jobs/SummaryStats.Mapping.v1.out
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -M edahi@lji.org
#PBS -l mem=1GB
#PBS -m ae
#PBS -W depend=afterok:\${CollectQsubs}

# >>>>>> Cat all mapping stats into the Analysis folder:

cat \$analysis/02.Mapping/\$name/Colnames.MappingStats.csv \$analysis/02.Mapping/*/MappingStats.csv > \$analysis/AllMappingStats.csv

echo Seems like this is going to work, baby!
EOT

qsub $Jobs/SummaryStats.Mapping.v1.sh

EOF

# >>> If selected, Run job
if [ "$run" = "1" ]; then
#  echo "Running Job"
 qsub $Jobs/BaseScript.WGSeq.BSdownload.v1.sh
fi

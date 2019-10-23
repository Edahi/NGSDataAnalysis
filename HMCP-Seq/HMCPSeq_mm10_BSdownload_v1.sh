#!/bin/bash

# HMCP-Seq BaseSpace Download Standarized pipeline mm10: Version -- 01
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2019.03.12
# Description: download fastq from Illumina's BaseSpace and analize data.
# Versions logs
# v1.
#  Download data, create scripts for each downloaded sample and run it
#  Current version linked to:
#  1) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPCseq_mm10_Existing_v1
#  2) 
#
# >>> How to use the program:
#    HMCPSeq_mm10_BSdownload_v1.sh [-b BaseSpace_project_ID || -t Run TroubleShooting] [-s Should someone do SE sequencing] [-d DownloadDir def:/BioScratch/edahi] [-a AnalysisDir def:/BioScratch/edahi] [-r Run created job def:no] [-q rao-exclusive queue def:no] [-h help]"
# >>> Example
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPSeq_mm10_BSdownload_v1.sh -b 122477370  -d /mnt/BioScratch/edahi/3_6_19_Jerry_Pool3 -a /mnt/BioAdHoc/Groups/RaoLab/Edahi/09.CombinedProjects/04.Jerry_Vipul_Hiroshi_Daniela/HMCP/3_6_19_Jerry_Pool3 -q -r
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPSeq_mm10_BSdownload_v1.sh -b 122471567  -d /mnt/BioScratch/edahi/3_6_19_Jerry_Pool4 -a /mnt/BioAdHoc/Groups/RaoLab/Edahi/09.CombinedProjects/04.Jerry_Vipul_Hiroshi_Daniela/HMCP/3_6_19_Jerry_Pool4 -q -r
# 

usage()
{
    printf "\nusage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPSeq_mm10_BSdownload_v1.sh "
    printf "\n\t[ -b | --BS        BaseSpace_project_ID | 'PROXI'   def:NONE -- Required Argument ]"
    printf "\n\t[ -d | --download  DownloadPath                     def:/mnt/beegfs ]"
    printf "\n\t[ -a | --analysis  AnalysisPath                     def:/mnt/BioScratch/edahi/ATAC ]"
    printf "\n\t[ -q | --queue     set 'rao-exclusive' queue        def:'default' ]"
    printf "\n\t[ -r | --run       Run created job                  def:no   ]"
    printf "\n\t[ -k | --keep      Keep intermediate results        def:no   ]"
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
analysis=/mnt/BioScratch/edahi/HMCP
     pid=
    name=HMCP
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

Jobs=$analysis/Jobs
mkdir -p $Jobs

# >>> Check if TroubleShooting:
# if [ "$ts" = "1" ]; then
#  cat <<EOF> $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.sh
# #!/bin/bash -ex
# #PBS -N TroubleShooting.HMCPSeq.BSdownload.mm10.v1
# #PBS -l walltime=168:00:00
# #PBS -o $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.out
# #PBS -e $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.out
# #PBS -j oe
# #PBS -l nodes=1:ppn=4
# #PBS -M edahi@lji.org
# #PBS -l mem=10GB
# #PBS -m ae
# EOF
# 
# if [ "$raoqueue" = "1" ]; then
#  cat <<EOF>> $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.sh
# #PBS -q rao-exclusive
# EOF
# 
# else
#  cat <<EOF>> $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.sh
# #PBS -q default
# EOF
# fi
# 
# cat <<EOF>> $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.sh
# From Script:
#       pid=$pid
#  analysis=$analysis
#  download=$download
#  raoqueue=$raoqueue
#      keep=$keep
#      Jobs=$Jobs
#      temp=$temp
# Obtain dataset's names:
# nameS=(\$(basename -a \$(ls $download/* | grep L00._R1 | rev | cut -c22- | rev | uniq)))
# Iterate throught the names of the downloaded sets AND
# Build the TroubleShooting script for each file:
# for i in \${!nameS[@]}; do 
#  name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g'); # echo \$name -- \${nameS[\$i]}
#  /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_Existing_TroubleShooting_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis -q -r $ts1 $ts2 $ts3
# done
# EOF
# 
#  if [ "$run" = "1" ]; then
#   echo "Running TroubleShooting For Whole Library"
#   qsub $Jobs/TroubleShooting.HMCPSeq.BSdownload.mm10.v1.sh
#  fi
#  exit
# fi

# >>>>>> Other than Troubleshooting:
# >>> Check if PID was given:
if [ "$pid" = "" ]; then
    printf "\n\tPID CODE is a required argument if not TroubleShooting. \nCheck --help for further assistance\n\n";
    exit
fi

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


# >>> Download BaseSpace data:
cat <<EOF> $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.sh
#!/bin/bash -ex
#PBS -N BaseScript.HMCPSeq.BSdownload.mm10.v1
#PBS -l walltime=168:00:00
#PBS -o $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.out
#PBS -e $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
EOF
if [ "$raoqueue" = "1" ]; then
 cat <<EOF>> $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.sh
#PBS -q rao-exclusive
EOF
else
 cat <<EOF>> $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.sh
#PBS -q default
EOF
fi
cat <<EOF>> $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.sh
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

 if [ "\$CollectQsubs" = "" ]; then
   CollectQsubs=\$(/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPSeq_mm10_Existing_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis $Parameters)

 else
   CollectQsubs=\${CollectQsubs}:\$(/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/HMCPSeq_mm10_Existing_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis $Parameters)
 fi
done

# >>> After all the jobs are completed:
# cat <<EOT> $Jobs/HMCPSeq.PostProcessing.v1.sh
# #!/bin/bash -x
# #PBS -N HMCPSeq-PostProcessing-v1
# #PBS -l walltime=8:00:00
# #PBS -o $Jobs/HMCPSeq.PostProcessing.v1.out
# #PBS -e $Jobs/HMCPSeq.PostProcessing.v1.out
# #PBS -j oe
# #PBS -l nodes=1:ppn=1
# #PBS -M edahi@lji.org
# #PBS -l mem=10GB
# #PBS -m ae
# #PBS -W depend=afterok:\${CollectQsubs}
# # >>>
# # >>> PLACE
# # >>> HERE
# # >>> NEXT
# # >>> STEPS
# # >>>
# 
# echo Seems like this is going to work, baby!
# EOT
# 
# qsub $Jobs/HMCPSeq.PostProcessing.v1.sh

EOF

# >>> If selected, Run job
if [ "$run" = "1" ]; then
 echo "Running Job"
 qsub $Jobs/BaseScript.HMCPSeq.BSdownload.mm10.v1.sh
fi

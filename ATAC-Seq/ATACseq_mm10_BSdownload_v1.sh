#!/bin/bash

# ATAC-seq BaseSpace Download Standarized pipeline mm10: Version -- 1.2
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2020.03.30
# Description: download fastq from Illumina's BaseSpace and analize data.
# Versions logs
# v1.2
#   Added PID to Job's name
#   Highlighted in red improper arguments
#   Added date to versioning
#   Added Date and Version to "Usage" printout
#   Updated HELP format
#   
# v1.1
#   Removed cluttered comment section at the beggining
#   Recorded the script name and arguments/parameters given for the log.
#   Updated line to check if arguments/parameters were given.
#   Updated cat's End Of File syntax for better readability.
#   Added pipeline versioning
#   Added version to generated scripts
#   Updated queue selection
#   Updated line to run generated shell scripts
#   Missing spaces in the usage function
#   Added "exit 1" as part of the usage function.
#   Moved up (after the parameters case selection) the parameters recollection for daughter jobs
#   Updated BS download perl script line
#   Removed unnecesary echo commands "Running jobs" and moved the trouble shooting one at the beggining of the if capsule
#   Updated prompted message if invalid argument was given
#   
# v1.0
#  Download data, create scripts for each downloaded sample and run it
#  Current version linked to:
#  1) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_Existing_v1
#  2) /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_Existing_TroubleShooting_v1

# >>> Pipeline Version:
 Ver=1.2
Date=2020.03.30

usage()
{
printf "
    Version:\t %s
       Date:\t %s
\n      Usage: %s    <PARAMETERS>\n
Req:\t[ -b | --BS            <#> ] BaseSpace projectID (Use 'PROXI' if data downloaded before)
\t[ -d | --download    <str> ]  Absolute <PATH> Folder to store downloaded data (def:/mnt/beegfs)
\t[ -a | --analysis    <str> ]  Absolute <PATH> Folder to conduct analyses (def:/mnt/BioAdHoc/Groups/RaoLab/temp)
\t[ -q | --queue             ]  Set 'rao-exclusive' queue   (def:'default')
\t[ -r | --run               ]  Run created job  (def:no)
\t[ -k | --keep              ]  Keep intermediate results        (def:no)
\t[ -t | --ts                ]  Run TroubleShooting for Library  (def:no)
\t[ -1 | --ts1               ]  Set TroubleShooting #1           (def:no)
\t[ -2 | --ts2               ]  Set TroubleShooting #2           (def:no)
\t[ -3 | --ts3               ]  Set TroubleShooting #3           (def:no)
\t[ -h | --help              ]  Show this message and exit                \n\n" $Ver $Date $0
exit 1
}

# >>> Check if parameters given:
[[ "$1" == "" ]] && usage # <<< Ver v1.1

# >>> Record parameters given: # <<< Ver v1.1
CodeLine=$0
CodeLine=${CodeLine}" "$@ 

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
                                ;;
        * )                     tput setaf 1 && tput bold && tput smul && printf "\n\tInvalid argument provided: %s\n\n" $1 && tput sgr0
		                        grep --color=always -e "$1" <<< "${CodeLine}" 
								usage
    esac
    shift
done

# >>> Obtain parameters for daughter script: # <<< Ver v1.1
[[ "$raoqueue" == "1" ]] && { raoqueue=rao-exclusive && Parameters=" "${Parameters}" -q" ; } || raoqueue=default
[[ "$run"      == "1" ]] && Parameters=${Parameters}" -r"
[[ "$keep"     == "1" ]] && Parameters=${Parameters}" -k"

# >>> Print Settings information:
printf '\n# >>> %s\nCommand:\n%s\n' $pid "${CodeLine}" # <<< Ver v1.2
printf '\nBaseSpace ID #:\t%s\nPBS Queue:\t%s\n\n' $pid $raoqueue
printf 'Parameters Config for daugther jobs: %s' "${Parameters[@]}" # <<< Ver v1.1
echo  ''

Jobs=$analysis/Jobs
mkdir -p $Jobs

# >>> Check if TroubleShooting:
if [ "$ts" = "1" ]; then
 printf "Running TroubleShooting\n" # <<< Ver v1.1
 cat > $Jobs/TroubleShooting.ATACseq.BSdownload.mm10.v${Ver}.sh <<EOF # <<< Ver v1.1
#!/bin/bash -ex
#PBS -N TroubleShooting.ATACseq.BSdownload.mm10.v${Ver}
#PBS -l walltime=24:00:00
#PBS -o $Jobs/TroubleShooting.ATACseq.BSdownload.mm10.v${Ver}.out
#PBS -e $Jobs/TroubleShooting.ATACseq.BSdownload.mm10.v${Ver}.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
#PBS -q ${raoqueue} # <<< Ver v1.1

# From Script:
      pid=$pid
 analysis=$analysis
 download=$download
     keep=$keep
     Jobs=$Jobs
     temp=$temp
# Obtain dataset's names:
nameS=(\$(basename -a \$(ls \$download/* | grep L001_R1)))
nameS=("\${nameS[@]/_L001_R1_001.fastq.gz/}")
# Iterate throught the names of the downloaded sets AND
# Build the TroubleShooting script for each file:
for i in \${!nameS[@]}; do 
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g');
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_Existing_TroubleShooting_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis -q -r $ts1 $ts2 $ts3
done
EOF

 [[ "$run" = "1" ]] && qsub $Jobs/TroubleShooting.ATACseq.BSdownload.mm10.v${Ver}.sh  # <<< Ver v1.1
 exit 1
fi

# >>>>>> Other than Troubleshooting:
# >>> Check if PID was given:
[[ "$pid" == "" ]] && { printf "\n\tProject ID (--BS) is a required argument if not TroubleShooting.\n\n" && usage ; } # Ver v1.1

# >>> Download BaseSpace data:
cat > $Jobs/BaseScript.ATACseq.BSdownload.mm10.v${Ver}.sh <<EOF # <<< Ver v1.1
#!/bin/bash -ex
#PBS -N BaseScript.ATACseq.BSdownload.mm10.v${Ver}.${pid}
#PBS -l walltime=24:00:00
#PBS -o $Jobs/BaseScript.ATACseq.BSdownload.mm10.v${Ver}.${pid}.out
#PBS -e $Jobs/BaseScript.ATACseq.BSdownload.mm10.v${Ver}.${pid}.out
#PBS -j oe
#PBS -l nodes=1:ppn=4
#PBS -M edahi@lji.org
#PBS -l mem=10GB
#PBS -m ae
#PBS -q ${raoqueue} # <<< Ver v1.1

# From Script:
      pid=$pid
 analysis=$analysis
 download=$download
     keep=$keep
     Jobs=$Jobs
     temp=$temp

# >>>>>> Download data from BaseSpace using MY CREDENTIALS (Project must be shared with me in order to work).
[[ "\$pid" = "PROXI" ]] && printf "Data downloaded before, skipping download...\n" || perl /share/apps/BSDownload/BSDownload-Generic.pl -c 739a72400375493aac50a8aeb910c3eb -s 7e4fd758cef94f31845a51116f007373 -t 7acae786b29b4a00ae883f6b80edaef5  -o \$download -p \$pid # <<< Ver v1.1

# Obtain dataset's names:
lane=\` ls \$download/*L00*_R*.gz | head -n1 | rev | cut -c17\`
nameS=(\$(basename -a \$(ls \$download/* | grep L00\${lane}_R1)))
nameS=("\${nameS[@]/_L00\${lane}_R1_001.fastq.gz/}")

# Iterate throught the names of the downloaded sets AND
# Build the analysis script for each file:
for i in \${!nameS[@]}; do 
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g');
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/ATACseq_mm10_Existing_v1 -s \${nameS[\$i]} -n \$name -d \$download -a \$analysis $Parameters
done
EOF

# >>> If selected, Run job
[[ "$run" = "1" ]] && qsub $Jobs/BaseScript.ATACseq.BSdownload.mm10.v${Ver}.sh

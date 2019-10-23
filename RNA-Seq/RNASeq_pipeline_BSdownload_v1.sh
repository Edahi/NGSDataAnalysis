#!/bin/bash

# RNA-Seq Super Pipeline for hg38|mm10 Single|Paired End data from one sequencing run: Version -- 01.5
# Author: Edahi Gonzalez-Avalos
# Email:  edahi@lji.org
# Date:   2019.09.30
# Versions logs
# v1.5
#   Rename tem files for reports to include the random seed to avoid overlap with other running programs
#   Print Version when help is displayed
#   Module'd the read length to 25-scale to never miss an STAR reference, unless read is bigger than 250...
#   Also, there is no longer a requirement to trim the reads, since they will get the closer reference possible
# v1.4
#   Reorder the "declare variable" section to match that of the argument reading
#   Inform user of the runQC and runDEG seceltions
#   Change "PBS queue" display information by "Parameters" variable content
# v1.3
#   Moved the folders declaration content inside the script for more generalizable content
#   Added a variable for the version
#   Renamed Jobs to include version used
# v1.2
#  Updated subread's featureCounts  subread-1.5.3-Linux-x86_64 ==> subread-2.0.0-source
# v1.1
#  Updated HOMER's functions used
#  Added optional Differential Expression analysis parameter
#  Re-Ordered the QC analysis and recolection to make independent FASTQC and conda TIN's
# v1.0
#  Mapping

# >>> Pipeline Version:
Ver=1.5

# >>> How to use the program:

usage()
{
printf "
\tVersion:\t%s\n
usage: /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/Bash/RNASeq_pipeline_BSdownload_v1.sh    <PARAMETERS>
\nReq:\t[-b | --BaseSpace   <#> ]  ProjectID (use 'PROXI' if downloaded before) def:NONE -- Required Argument 
\nReq:\t[-g | --Genome    <str> ]  Genome Version  mm9|mm10|hg18|hg19|hg38      def:NONE -- Required Argument 
\nReq:\t[-c | --Config    <str> ]  Differential Expression Configurations       def:NONE -- Required Argument 
\n\t[-d | --download  <str> ]  DownloadPath                                 def:/mnt/beegfs 
\n\t[-a | --analysis  <str> ]  AnalysisPath                                 def:/mnt/BioScratch/edahi/ 
\n\t[-m | --Mapper    <str> ]  Mapping Tool    STAR|tophat2                 def:STAR 
\n\t[-s | --MapStrand   <#> ]  0:Unstanded|1:Stranded|2:ReverselyStranded   def:0 
\n\t[-t | --temp      <str> ]  TemporalPath                                 def:/mnt/BioScratch/edahi 
\n\t[-b | --barcode     <#> ]  Barcode sequence length to remove --if any   def:0 
\n\t[-q | --queue           ]  Set 'rao-exclusive' queue                    def:'default' 
\n\t[-r | --run             ]  Run generated job                            def:no   
\n\t[-k | --keep            ]  Keep intermediate results                    def:no   
\n\t[     --noQC            ]  Turn Off QC analysis                         def:On
\n\t[     --noDEG           ]  Turn Off Differential Expression analysis    def:On
\n\t[     --NoAnal          ]  Stop after Data download                     def:Off
\n\t[-h | --help            ]  Show this message and exit\n\n" $Ver
}

# >>> Check if arguments given:
if [ "$1" == "" ]; then
  usage
  exit 1
fi

# >>> Declare Variables
download=/mnt/beegfs
analysis=/mnt/BioScratch/edahi
Mapper=STAR
MapStrand=0
surname=
raoqueue=0
run=0
keep=0
runQC=1
DEGanal=1
barcode=0
TEMP=/mnt/BioScratch/edahi
RandomSeed=$RANDOM

# >>> Assign arguments to variables
while [ "$1" != "" ]; do
    case $1 in
        -b | --BS )             shift
                                pid=$1
                                ;;
        -g | --GenVer )         shift
                                GenVer=$1
                                ;;
        -c | --Config )         shift
                                Config=$1
                                ;;
        -d | --download )       shift
                                download=$1
                                ;;
        -a | --analysis )       shift
                                analysis=$1
                                ;;
        -m | --Mapper )         shift
                                Mapper=$1
                                ;;
        -s | --MapStrand  )     shift
                                MapStrand=$1
                                ;;
        -q | --queue  )         raoqueue=1
                                ;;
        -r | --run    )         run=1
                                ;;
        -k | --keep  )          keep=1
                                ;;
             --noQC  )          runQC=0
                                ;;
             --noDEG )          DEGanal=0
                                ;;
        -b | --barcode )        shift
                                barcode=$1
                                ;;
        -t | --temp )           shift
                                temp=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     printf "\n\tAn invalid argument was given: %s \n" $1
                                usage
                                exit 1
    esac
    shift
done

# >>> Check if File surname was given:
[[ $pid = "" ]] &&  printf "\n\tBasespace's Project ID '--BS' is a required argument if not TroubleShooting. \n\tCheck --help for further assistance.\n\n" && exit

# >>> Check if valid Genome version build was given:
list=(mm9 mm10 hg18 hg19 hg38)
Over(){ printf "\n\tProvide a valid genome build version \n\tCheck --help for further assistance\n\n" ; exit;}
[[ ${list[@]} =~ (^|[[:space:]])$GenVer($|[[:space:]]) ]]  ||  Over
[[ $GenVer =~ (^hg) ]] && Organism=Homo_sapiens ||  Organism=Mus_musculus

# >>> Check if valid Aligner was given:
list=(STAR tophat2)
Over(){ printf "\n\tProvide a valid RNA-Seq Aligner \n\tCheck --help for further assistance\n\n" ; exit;}
[[ ${list[@]} =~ (^|[[:space:]])$Mapper($|[[:space:]]) ]]  ||  Over

# >>> Check if valid Strandness was given:
list=(0 1 2)
Over(){ printf "\n\tProvide a valid RNA-Seq Strandness \n\tCheck --help for further assistance\n\n" ; exit;}
[[ ${list[@]} =~ (^|[[:space:]])$MapStrand($|[[:space:]]) ]]  ||  Over

# >>> Obtain parameters for script:
Parameters=
[[ $raoqueue == "1" ]]  && Parameters=" "${Parameters}" -q" && raoqueue=rao-exclusive ||  raoqueue=default
[[ $run == "1" ]]       && Parameters=" "${Parameters}" -r"
[[ $keep == "1" ]]      && Parameters=" "${Parameters}" -k"

# >>> Print Settings information:
printf '\nBS ProjectID:\t%s\nParameters:\t%s\nOrganism:\t%s \nGenome Version:\t%s\nMappingTool:\t%s\nStrandness:\t%s\nRandomSeed:\t%s\nRun QC:\t\t%s\nRun DiffExp:\t%s\n\n' $pid $Parameters $Organism $GenVer $Mapper $MapStrand $RandomSeed $runQC $DEGanal

Jobs=$analysis/Jobs
mkdir -p $Jobs


if [ "$DEGanal" = "1" ]; then

# >>> Generating Script for Differential Gene expression:
# cat >>  $Jobs/CMSIP.${GenVer}.SRR.v1.${name}.sh << EOF

cat > $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R <<EOF
#!/share/apps/R/3.5.2/bin/Rscript

library(IRdisplay)
library(limma)
library(edgeR)
library(Glimma)
library(Mus.musculus)
library(RColorBrewer)
library(gplots)
library(clusterProfiler)
library(pathview)
library(DOSE)
library(TCseq)


setwd("$analysis/06.DiffExp")
counts <- read.table("$analysis/04.FeatureCounts/${Mapper}_gene_counts_s${MapStrand}.txt", header=TRUE)
rownames(counts) <- as.vector(counts[,1])
counts <- counts[,6:ncol(counts)]

# >>> Name Condition(s):
Configs     = read.table("$Config", header=F, stringsAsFactors=FALSE)
condition   = c(Configs[,1])
UniqueNames = colnames(counts)[-1]
UniqConfigs = unique(Configs[,1])

# >>> Use main Condition to generate the DGE object
DGE = DGEList(counts[,2:(length(condition)+1)], group = condition)
samplenames = colnames(DGE\$counts)
DGE\$samples\$UniqueSamples = as.factor(UniqueNames)

# >>> Obtain EntrezID for future annotation and add to main dataset:
geneid = rownames(DGE)
genes  = select(Mus.musculus, keys=geneid, columns=c("ENTREZID", "TXCHROM"), keytype="SYMBOL")
# genes  = genes[!duplicated(genes\$ENTREZID),]
genes  = genes[match(geneid,genes[,1]),]
DGE\$genes = genes

# >>> Normalization:
# For differential expression and related analyses, gene expression is rarely considered at the level of raw counts since libraries sequenced at a greater depth will result in higher counts. 
#   Rather, it is common practice to transform raw counts onto a scale that accounts for such library size differences. 
#   Popular transformations include counts per million (CPM), log2-counts per million (log-CPM), reads per kilobase of transcript per million (RPKM), and fragments per kilobase of transcript per million (FPKM).
# In our analyses, CPM and log-CPM transformations are used regularly although they do not account for feature length differences which RPKM and FPKM values do. 
#   Whilst RPKM and FPKM values can just as well be used, CPM and log-CPM values can be calculated using a counts matrix alone and will suffice for the type of comparisons we are interested in. 
#   Assuming that there are no differences in isoform usage between conditions, differential expression analyses look at gene expression changes between conditions rather than comparing expression across multiple genes or drawing conclusions on absolute levels of expression. 
#   In other words, gene lengths remain constant for comparisons of interest and any observed differences are a result of changes in condition rather than changes in gene length.
# Here, raw counts are converted to CPM and log-CPM values using the cpm function in edgeR, where log-transformations use a prior count of 0.25 to avoid taking the log of zero. 
#   RPKM values are just as easily calculated as CPM values using the rpkm function in edgeR if gene lengths are available.

# > Perform the CPM Normalization as backup for future plot showing differences Before-After Normalization.
cpm  = cpm(DGE)
lcpm = cpm(DGE, log=TRUE)
# > See How many genes does not have a single raw count in any of the comparisons
table(rowSums(DGE\$counts==0)==(dim(DGE\$counts)[2]))
# > index for genes whose at least a third out of my samples have 1 count per million
keep.exprs = rowSums(cpm>1)>=round((dim(DGE\$counts)[2])/3) 

# The lower the coverage, the stringent the threshold 
# "Because you wont trust if a lowly expressed gene is actually being expressed or not"
DGE = DGE[keep.exprs,, keep.lib.sizes=FALSE]
dim(DGE)

# This cell supports future resizing of the plots in the Jupyter notebook
defaultWidth  = getOption("repr.plot.width")
defaultHeight = getOption("repr.plot.height")
# This cell resizes the output for attractive side-by-side plots
options(repr.plot.width=defaultWidth, repr.plot.height=defaultHeight*.7)

# >>> Display Normalization changes:
ColorSets=c(brewer.pal(name="Paired", n = 12),brewer.pal(name="Dark2", n = 8))
nsamples = ncol(DGE)
col      = ColorSets[1:nsamples]
pdf("Count_Data_Normalization_Graphs_1_Count_Density_Plots.pdf")
 par(mfrow=c(1,2))
 plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
  title(main="Normalized Unfiltered Count Data", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){ den <- density(lcpm[,i]) ;  lines(den\$x, den\$y, col=col[i], lwd=2)}
  legend("topright", samplenames, text.col=col, bty="n")
  lcpm = cpm(DGE, log=TRUE)
 plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
  title(main="Normalized Filtered Count Data", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){ den <- density(lcpm[,i]) ;   lines(den\$x, den\$y, col=col[i], lwd=2)}
  legend("topright", samplenames, text.col=col, bty="n")
dev.off()


# >>> Normalization proof:
DGE    = calcNormFactors(DGE, method = "TMM")
DGE\$samples\$norm.factors
DGE.2  = DGE
DGE.2\$samples\$norm.factors = 1
DGE.2\$counts[,1] = ceiling(DGE.2\$counts[,1]*0.05)
DGE.2\$counts[,2] = DGE.2\$counts[,2]*5
pdf("Count_Data_Normalization_Graphs_2_Count_BoxPlots.pdf")
 lcpm = cpm(DGE.2, log=TRUE)
 boxplot(lcpm, las=2, col=col, main="",ylim=c(-6,17))
  title(main="Unnormalized data",ylab="Log-cpm")
 DGE.2 = calcNormFactors(DGE.2)  
 DGE.2\$samples\$norm.factors
 lcpm = cpm(DGE.2, log=TRUE)
 boxplot(lcpm, las=2, col=col, main="")
  title(main="Normalized data",ylab="Log-cpm")
dev.off()

# >>> Samples Clustering:
lcpm = cpm(DGE, log=TRUE)
col.group <- as.factor(condition)
col.sample <- as.factor(UniqueNames )
levels(col.group) =    ColorSets[1:nlevels(col.group)]
levels(col.sample)   =    ColorSets[1:nlevels(col.sample)]
col.group    = as.character(col.group)
col.sample   = as.character(col.sample)
pdf("Count_Data_Normalization_Graphs_3_MultiDimensionalScaling.pdf")
 plotMDS(lcpm, labels=condition, col=col.group)
  title(main="A. Sample groups")
 plotMDS(lcpm, labels=UniqueNames, col=col.sample)
  title(main="B. Individual Samples")
dev.off()
glMDSPlot(lcpm, labels=colnames(counts)[-1], folder = "MultiDimensionalScaling_Explorative", groups=DGE\$samples[,c(1,4)])

# >>> Differential Analysis Design:
Treat  = DGE\$samples\$group
design = model.matrix(~0+Treat)
colnames(design) = levels(Treat)
fit  = lmFit(DGE\$counts,design)
# Now we can make any comparisons between the experimental conditions:
contr.matrix <- makeContrasts(
EOF

UniqCond=(`cut -f1 $Config  | sort | uniq `)
for i in ${!UniqCond[@]}; do 
 for j in `seq ${i} $((${#UniqCond[@]}-1))`; do 
  [[ $i != $j ]] && echo ${UniqCond[${i}]}_${UniqCond[${j}]} = ${UniqCond[${i}]}-${UniqCond[${j}]}, >> $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R
 done
done

cat >> $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R <<EOF
 levels=design)

contr.matrix
pdf("Count_Data_Normalization_Graphs_3_Variance.pdf")
 par(mfrow=c(1,2))
 v = voom(DGE, design, plot=TRUE)
 v
 options(repr.plot.width=defaultWidth, repr.plot.height=defaultHeight)
 vfit <- lmFit(v, design)
 vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
 efit <- eBayes(vfit)
 plotSA(efit, main="Final model: Mean-variance trend")
dev.off()

# For a quick look at differential expression levels, the number of significantly up- and down-regulated genes can be summarised in a table. 
#   Significance is defined using an adjusted p-value cutoff that is set at 5% by default. 
#   1st exp) For the comparison between expression levels in StimuliUS vs StimuliPI for Plus, 745  genes are found to be down-regulated in US relative to PI and 724  genes are up-regulated in US vs PI for Plus â€“ a total of 1,469 DE genes. 
#   The larger numbers of DE genes observed for comparisons involving the basal population are consistent with our observations from the MDS plots.
summary(decideTests(efit))

#1st Experiment: A_vs_B_for_x   -1 = Lower in A   +1 = Higher in A  >|<  
#  # PlusvsMinusForStimuliUS PlusvsMinusForStimuliPI StimuliUSvsStimuliPIforMinus   StimuliUSvsStimuliPIforPlus
# -1                      18                      25                          681                           745
# 0                    10941                   10953                         9551                          9604
# 1                      114                      95                          841                           724

# >>> Apply a log fold change threshold to the Differential Analysis:
tfit = treat(vfit, lfc=log2(1))
dt   = decideTests(tfit)
summary(dt)

# >>> Paired scatterplots comparisons:
pairs.raw  = counts[keep.exprs,-1]
pairs.lcpm = lcpm
pairs.v.e  = v\$E        # <<< numeric matrix of normalized expression values on the log2 scale
pairs.v.w  = v\$weights  # <<< numeric matrix of inverse variance weights

# >>>Spearman's
panel.tab <- function (x,y) {
  par(new = TRUE);
  smoothScatter(x,y); 
  legend(legend=c(paste0("cor=",round(cor(x,y,method="spearman"),2))),"topleft",cex=1.25,bty="n")}
pdf("Count_Data_Normalization_Graphs_4_Spearman_Correlations.pdf",width=14,height=14)
 pairs(pairs.raw , panel = panel.tab, lower.panel = NULL, main = "Raw Counts")
 pairs(pairs.lcpm, panel = panel.tab, lower.panel = NULL, main = "Log2 Counts Per Million")
 pairs(pairs.v.e , panel = panel.tab, lower.panel = NULL, main = "Voom Normalized expression values on the log2 scale")
 pairs(pairs.v.w , panel = panel.tab, lower.panel = NULL, main = "Voom Inverse Variance Weights")
dev.off()

# >>>Pearson's
panel.tab <- function (x,y) {
  par(new = TRUE);
  smoothScatter(x,y); 
  legend(legend=c(paste0("cor=",round(cor(x,y,method="pearson"),2))),"topleft",cex=1.25,bty="n")}
pdf("Count_Data_Normalization_Graphs_4_Pearson_Correlations.pdf",width=14,height=14)
 pairs(pairs.raw , panel = panel.tab, lower.panel = NULL, main = "Raw Counts")
 pairs(pairs.lcpm, panel = panel.tab, lower.panel = NULL, main = "Log2 Counts Per Million")
 pairs(pairs.v.e , panel = panel.tab, lower.panel = NULL, main = "Voom Normalized expression values on the log2 scale")
 pairs(pairs.v.w , panel = panel.tab, lower.panel = NULL, main = "Voom Inverse Variance Weights")
dev.off()

# de.common = which(dt[,1]!=0 & dt[,2]!=0)
# length(de.common)
# head(tfit$genes$SYMBOL[de.common], n=20)
# pdf("DiffExpression_Graphs_1_Visual_Representation_of_Shared_Hits.pdf")
# 
# # >>> correct
#  vennDiagram(dt[,1:3], circle.col=c("turquoise", "salmon", "Green"),names=c("Rest_CsA","Rest_CsAIono","Rest_Iono"), main="\n\nDifferentially Expressed Genes\nOverlapped Across Comparisons.")
# dev.off()
# write.fit(tfit, dt, file="DEGenes_adjPval_0.05_Log2FC_1.txt")
# # >>> topTread:      Extract a table of the top-ranked genes from a linear model fit.



# >>> Heatmaps of the most Variable Genes
mycol <- colorpanel(1000,"blue","white","red")

EOF

coef=0
for i in ${!UniqCond[@]}; do 
 for j in `seq ${i} $((${#UniqCond[@]}-1))`; do 
  [[ $i != $j ]] && coef=$(($coef+1))
  [[ $i != $j ]] && echo '# >>> Interactive Differential Expression Plots
'${UniqCond[${i}]}_${UniqCond[${j}]} = topTreat\(tfit, coef=$coef, n=Inf\)                                           '
'pdf\(\"DiffExpression_Graphs_2_MAPlot_${UniqCond[${i}]}-${UniqCond[${j}]}.pdf\"\)                                   '     
' plotMD\(tfit, column=$coef, status=dt[,$coef], main=colnames\(tfit\)[$coef], xlim=c\(-8,13\), ylim=c\(-15,15\)\)   '
'dev.off\(\)                                                                                                         '
'glMDPlot\(tfit, coef=$coef, status=dt, main=colnames\(tfit\)[$coef],side.main=\"SYMBOL\", counts=DGE\$counts,       '
' groups=condition,  folder=\"Explorative_Differential_Expression_${UniqCond[${i}]}-${UniqCond[${j}]}\"\)            '

# >>> Heatmap of most variable Genes
'${UniqCond[${i}]}_${UniqCond[${j}]}.top      = ${UniqCond[${i}]}_${UniqCond[${j}]}\$SYMBOL[1:100]                   '
'i = which\(v\$genes\$SYMBOL %in% ${UniqCond[${i}]}_${UniqCond[${j}]}.top\)[1:100]                                   '
'pdf\(\"DiffExpression_Graphs_3_Heatmap_100_${UniqCond[${i}]}-${UniqCond[${j}]}.pdf\"\)                              '
' options\(repr.plot.width=defaultWidth, repr.plot.height=defaultHeight*1.5\)                                        '
' heatmap.2\(v\$E[i,], scale=\"row\", labRow=v\$genes\$SYMBOL[i][1:100], labCol=condition,                           '
'  col=mycol, trace=\"none\", density.info=\"none\",                                                                 '
'  margin=c\(8,6\), lhei=c\(2,10\), dendrogram=\"column\",cexRow=0.4,main = \"${UniqCond[${i}]}-${UniqCond[${j}]}\"\)'
'dev.off\(\)                                                                                                         '

# >>> Genome Onthology 
'ToGo=cbind\(${UniqCond[${i}]}_${UniqCond[${j}]}\$logFC , ${UniqCond[${i}]}_${UniqCond[${j}]}\$adj.P.Val\)           '
'colnames\(ToGo\)=c\(\"logFC\",\"adjPVal\"\)                                                                         '
'rownames\(ToGo\)=${UniqCond[${i}]}_${UniqCond[${j}]}\$ENTREZID                                                      '
'ToGo = ToGo[!is.na\(rownames\(ToGo\)\),]                                                                                '
'gene=rownames\(ToGo\)[ToGo[,1] \>  2]                                                                                '
# >>> GO over-representation Upregulated
'${UniqCond[${i}]}_${UniqCond[${j}]}_EGO_Up=enrichGO\(gene=gene, universe=rownames\(ToGo\), OrgDb=org.Mm.eg.db,      '
' ont=\"MF\", pAdjustMethod=\"BH\", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE\)                            '
'head\(${UniqCond[${i}]}_${UniqCond[${j}]}_EGO_Up\)                                                                  '
'gene=rownames\(ToGo\)[ToGo[,1] \< -2]                                                                                '
'gene=gene[!is.na\(gene\)]                                                                                           '
# >>> GO over-representation Downregulated
'${UniqCond[${i}]}_${UniqCond[${j}]}_EGO_Down=enrichGO\(gene=gene, universe=rownames\(ToGo\), OrgDb=org.Mm.eg.db,    '
' ont=\"MF\", pAdjustMethod=\"BH\", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE\)                            '
'head\(${UniqCond[${i}]}_${UniqCond[${j}]}_EGO_Down\)                                                                '


' >> $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R
 done
done
echo 'save(list = ls(all.names = TRUE), file = "'$Jobs/RNAseq.${GenVer}.${Mapper}.s${MapStrand}.v1'.RData", envir = .GlobalEnv)' >> $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R
echo '
# https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf
# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
# https://bioconductor.org/packages/release/bioc/html/DOSE.html
# https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf' >> $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R

fi #<<<<<< If DEG ( DEGanal ) is activated

# >>> Download BaseSpace data:
cat > $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.sh <<EOF
#!/bin/bash -ex
#PBS -N RNAseq.${GenVer}.${Mapper}.s${MapStrand}.v1
#PBS -l walltime=72:00:00
#PBS -o $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.out
#PBS -e $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.out
#PBS -j oe
#PBS -l nodes=1:ppn=8
#PBS -M edahi@lji.org
#PBS -l mem=20GB
#PBS -m ae
#PBS -q $raoqueue

# Prepare for QC analysis:
     runQC=$runQC


# From Script:
       pid=$pid
   barcode=$barcode
    GenVer=$GenVer
  Organism=$Organism
    Mapper=$Mapper
 MapStrand=$MapStrand
      TEMP=${TEMP}/RNA_analysis_rand${RandomSeed}
  download=$download
  analysis=$analysis
      Jobs=\$analysis/Jobs
  softlink=\$analysis/01.Data
       Map=\$analysis/02.Mapping
   TagDirs=\$analysis/03.TagDirectories
     FeatC=\$analysis/04.FeatureCounts
   BigWigs=\$analysis/05.BigWigs
       DEG=\$analysis/06.DiffExp
        QC=\$analysis/07.RSeQC
     runQC=$runQC
   DEGanal=$DEGanal
    Tracks=$analysis/Tracks
    FASTQC=$analysis/FASTQC
   multiQC=/home/edahi/download/code/MultiQC/scripts/multiqc
   
   
#Variables & files:
    fastqc=/Bioinformatics/apps/FastQC/FastQC-0.11.2/fastqc
    tophat=/home/danielasc/software/tophat-2.1.1.Linux_x86_64/tophat2
       bwa=/home/edahi/download/code/bwa/bwa
      STAR=/home/edahi/download/code/STAR/bin/Linux_x86_64_static/STAR
 mMultiWig=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/HOMER/bin/makeMultiWigHub.pl
   mTagDir=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/HOMER/bin/makeTagDirectory
    fcount=/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/subread-2.0.0-source/bin/featureCounts
         R=/share/apps/R/3.5.2/bin/R
       sam=/usr/bin/samtools

# >>>>>> Build directories:
mkdir -p \$FASTQC \$softlink \$Map \$TagDirs \$FeatC \$BigWigs \$Tracks \$TEMP \$DEG 

# >>>>>> Download data from BaseSpace using MY CREDENTIALS (Project must be shared with me in order to work).
[[ \$pid == "PROXI" ]] && printf "\nData downloaded before, skipping download...\n" ||  perl /share/apps/BSDownload/BSDownload-Generic.pl -c 739a72400375493aac50a8aeb910c3eb -s 7e4fd758cef94f31845a51116f007373 -t 7acae786b29b4a00ae883f6b80edaef5  -o \$download -p \$pid

# >>> Obtain dataset's names:
lane=\` ls \$download/*L00*_R*.gz | head -n1 | rev | cut -c10\`
nameS=(\$(basename -a \$(ls \$download/* | grep L00\${lane}_R1)))
nameS=("\${nameS[@]/_L00\${lane}_R1_001.fastq.gz/}")

for i in \${!nameS[@]}; do 
# >>> Get sample's given name:
 name=\$(echo \${nameS[\$i]} | sed 's@_S[0-9]\\+\$@@g');
 printf "\nSampleName:\t%s\n" \$name

# >>> check if PE:
 test -f \$download/\${nameS[\${i}]}_L00\${lane}_R2_*gz && PE=1 || PE=0
[[ \$PE == "1" ]] && layout=PE || layout=SE

# >>> Read length:
ReadLength=\$(zcat \$download/\${nameS[\$i]}_L00\${lane}_R1*fastq.gz --quiet | head -n2 | tail -n1 | tr -d '\n' | wc -c)
ReadLength=\`echo \$((\$ReadLength - \$barcode))\`

# >>> Determine STAR reference to use:
[[ \$PE == "1" ]] && STARref=\`echo \$((\$ReadLength * 2 ))\` || STARref=\$ReadLength
   diff=\`echo \$((\$STARref % 25))\`
if [ "\$diff" -lt 13 ]; then                    # <<< "Closer 25th STAR reference is Lesser than 13 units"
 echo \$((\$STARref - \$diff))
 STARref=\`echo \$((\$STARref - \$diff))\`
else                                            # <<< "Closer 25th STAR reference is Higher than 13 units"
 echo \$((\$STARref - \$diff + 25))
 STARref=\`echo \$((\$STARref - \$diff + 25))\`
fi

# >>> Merge lanes & rm barcodes (if any):
 if [ "\$barcode" = "0" ]; then
                        zcat \$download/\${nameS[\$i]}_L00*_R1_*fastq.gz  > \$TEMP/\${name}_R1.fastq && ln -s \$TEMP/\${name}_R1.fastq \$softlink/\${name}_R1.fastq
   [[ \$PE == "1" ]]  && zcat \$download/\${nameS[\$i]}_L00*_R2_*fastq.gz  > \$TEMP/\${name}_R2.fastq && ln -s \$TEMP/\${name}_R2.fastq \$softlink/\${name}_R2.fastq
 else
   printf "Removing %s Barcode bases...\n" \$barcode
                        zcat \$download/\${nameS[\$i]}_L00*_R1_*fastq.gz  | sed '2~2s/^.\{'\$barcode'\}//g'  > \$TEMP/\${name}_R1.fastq && ln -s \$TEMP/\${name}_R1.fastq \$softlink/\${name}_R1.fastq
   [[ \$PE == "1" ]]  && zcat \$download/\${nameS[\$i]}_L00*_R2_*fastq.gz  | sed '2~2s/^.\{'\$barcode'\}//g'  > \$TEMP/\${name}_R2.fastq && ln -s \$TEMP/\${name}_R2.fastq \$softlink/\${name}_R2.fastq
 fi

# >>> FASTQC
   \$fastqc \$softlink/\${name}_R1.fastq --outdir=\$FASTQC &
   [[ \$PE == "1" ]]  &&  \$fastqc \$softlink/\${name}_R2.fastq --outdir=\$FASTQC &

# >>> Map data  -- STAR --OR-- tophat2
           Genes=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Annotation/Archives/archive-current/Genes/genes.gtf
         starDir=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Sequence/STAR_ref_\${STARref}
      BowTie2Ref=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Sequence/Bowtie2Index/genome
 TranscriptIndex=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Sequence/TopHat2/TopHat2_transcript_index
          RefSeq=/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Annotation/Genes/RefSeq.union.bed
 [[ \$PE == "1" ]] && R2=\$softlink/\${name}_R2.fastq || R2= 
 [[ \$PE == "1" ]] && mMultiPar=" -flip -sspe"  || mMultiPar=" "

# >>> Map to tRNA / rRNA Followed by # >>> Index tRNA and rRNA mapping results
\$bwa mem -M -t 8 /mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Sequence/AbundantSequences/bwa/tRNA_\${GenVer}.fa \$softlink/\${name}_R1.fastq \$R2  > \$TEMP/\${name}.tRNA.sam && \
 \$sam view -uSh -F 4 \$TEMP/\${name}.tRNA.sam | \$sam sort -@ 8 - \$TEMP/\${name}.tRNA && \
 \$sam index     \$TEMP/\${name}.tRNA.bam && \
 \$sam flagstat  \$TEMP/\${name}.tRNA.bam > \$TEMP/\${name}.tRNA.Summary.txt &
\$bwa mem -M -t 8 /mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/\$Organism/UCSC/\$GenVer/Sequence/AbundantSequences/bwa/rRNA_\${GenVer}.fa \$softlink/\${name}_R1.fastq \$R2  > \$TEMP/\${name}.rRNA.sam && \
 \$sam view -uSh -F 4 \$TEMP/\${name}.rRNA.sam | \$sam sort -@ 8 - \$TEMP/\${name}.rRNA && \
 \$sam index     \$TEMP/\${name}.rRNA.bam && \
 \$sam flagstat  \$TEMP/\${name}.rRNA.bam > \$TEMP/\${name}.rRNA.Summary.txt &



# >>> Map to transcripts
 [[ \$Mapper == "STAR" ]]  &&         \$STAR      \\
  --genomeDir                     \$starDir      \\
  --outFileNamePrefix         \$Map/\${name}      \\
  --outFilterMultimapNmax               30      \\
  --outReadsUnmapped                 Fastx      \\
  --outSAMattributes                   All      \\
  --outSAMprimaryFlag         OneBestScore      \\
  --outSAMstrandField          intronMotif      \\
  --outSAMtype      BAM SortedByCoordinate      \\
  --quantMode                   GeneCounts      \\
  --readFilesIn  \$softlink/\${name}_R1.fastq \$R2 \\
  --readFilesCommand                   cat      \\
  --runRNGseed                        1921      \\
  --runThreadN                           8      \\
  --sjdbOverhang      \$((\$STARref -1))   ||     \\
                                                \\
  \$tophat                                       \\
  --transcriptome-index \$TranscriptIndex        \\
  -o                        \$Map/\${name}        \\
  -N                                   2        \\
  --read-gap-length                    1        \\
  --max-multihits                      1        \\
  -p                                   4        \\
  --no-coverage-search                          \\
  \$BowTie2Ref                                   \\
  \$softlink/\${name}_R1.fastq \$R2

# >>> Create Tag Directory: -- STAR --OR-- tophat2
#    -flip -sspe             For Paired-End Data
#    I should remove -sspe   For single-End Data
[[ \$Mapper == "STAR" ]]  && MapResults="Aligned.sortedByCoord.out.bam" || MapResults="/accepted_hits.bam"

\$mTagDir \$TagDirs/\${name}_\${Mapper} \$Map/\${name}\${MapResults} \\
    -fragLength given \\
    -format sam       \\
    -unique           \\
    -keepOne          \\
    -genome \$GenVer   \\
    \$mMultiPar        \\
    -checkGC 


# >>> MultiQC
cd \$FASTQC
\$multiQC . &

# >>>>>> Run QC analyses:
if [ "\$runQC" = "1" ]; then

 source /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/etc/profile.d/conda.sh
 conda activate RSeQC_py3
 mkdir -p \$QC/inner_distance \$QC/insertion_profile \$QC/junction_annotation \$QC/junction_saturation \$QC/infer_experiment \$QC/mismatch_profile \$QC/geneBody_coverage \$QC/Transcript_Integrity_Number

 # >>> Recollect Mapping results in a file:
 echo \$Map/\${name}\${MapResults} >> \$Map/MappingNames.txt
 
 # >>> Index Mapping Results:
 \$sam index  \$Map/\${name}\${MapResults}

 # >>> Calculate inner_distance
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/inner_distance.py      -i \$Map/\${name}\${MapResults} -r \$RefSeq -u 500 -o \$QC/inner_distance/\${name}       &   # For

 # >>> Calculate insertion_profile
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/insertion_profile.py   -i \$Map/\${name}\${MapResults} -s \$layout        -o \$QC/insertion_profile/\${name}    &   # For

 # >>> Calculate junction_annotation
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/junction_annotation.py -i \$Map/\${name}\${MapResults} -r \$RefSeq        -o \$QC/junction_annotation/\${name}  &   # For

 # >>> Calculate junction_saturation
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/junction_saturation.py -i \$Map/\${name}\${MapResults} -r \$RefSeq        -o \$QC/junction_saturation/\${name}  &   # For

 # >>> infer_experiment
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/infer_experiment.py    -i \$Map/\${name}\${MapResults} -r \$RefSeq         > \$QC/infer_experiment/\${name}     &   # For

 # >>> Calculate mismatch_profile
 /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/mismatch_profile.py    -i \$Map/\${name}\${MapResults} -l \$ReadLength    -o \$QC/mismatch_profile/\${name}     &   # For
fi

done

# ! # >>> If STAR Mapper, remove from shared memory Shared genome did not work
#         [[ \$Mapper == "STAR" ]]  &&  \$STAR --genomeDir \$starDir  --genomeLoad Remove

# >>> Generate Multi Wig Tracks:
\$mMultiWig RNAseq_\${Mapper} \${GenVer} \
 -url http://informaticsdata.liai.org/NGS_analyses/ad_hoc/\${Tracks/\/mnt\/BioAdHoc\//}/ \
 -webdir \$Tracks \
 -d \`ls \$TagDirs/* -d\`

# >>> Link normalized BigWig files:
cp \$Tracks/RNAseq_\${Mapper}/\${GenVer}/*bigWig \$BigWigs/ &

# >>> FeatureCounts. -- STAR --OR-- tophat2
cd \$FeatC/
[[ \$Mapper == "STAR" ]]  && MapResults="*Aligned.sortedByCoord.out.bam" || MapResults="*/accepted_hits.bam"

\$fcount   -o \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.txt \\
  -a \$Genes        \\
  -g gene_id       \\
  -s \$MapStrand    \\
  -p -B -T 8       \\
  \`ls \$Map/\${MapResults}\` > \$FeatC/\${Mapper}_s\${MapStrand}_Log.final.out 2>&1
             

# >>>>>> Mapping stats recollection:

# >>> Header:
echo "SampleName,TotalReads,MappedReads,AssignedReads_"\$Mapper"_s"\$MapStrand,tRNA_Reads,rRNA_Reads > \$FeatC/ReadAssignmentReport_HEAD_\${Mapper}_s\${MapStrand}.csv

# >>> Sample Names:
echo \${nameS[@]} | tr ' ' '\n' > \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv

# >>> Initial Reads:
[[ \$Mapper == "STAR" ]] && grep "Number of input" \$Map/*Log.final.out | rev | cut -f1 | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
if [ "\$Mapper" = "tophat2" ]; then
 if [ "\$PE" = "1" ]; then
  grep -e "Left reads" \$Map/*/align_summary.txt -A1 | grep Input | rev | cut -f1 -d\  | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
 else
  grep -e "Input" \$Map/*/align_summary.txt | rev | cut -f1 -d\  | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
 fi
fi

# >>> Mapped Reads
[[ \$Mapper == "STAR" ]] && grep -e "Uniquely mapped reads number" \$Map/*Log.final.out | rev | cut -f1  | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
if [ "\$Mapper" = "tophat2" ]; then
 if [ "\$PE" = "1" ]; then
  grep Aligned \$Map/*/align_summary.txt  | rev | cut -f1 -d\   | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
 else
  grep -e "Mapped" \$Map/*/align_summary.txt | cut -f1 -d\(  | rev | cut -f2 -d\   | rev | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
 fi
fi

# >>> Assigned Reads throught Feature Counts:
cat \$FeatC/\${Mapper}_s\${MapStrand}_Log.final.out   | grep -e "Success" | cut -f2 -d: | tr ' ' '\t' | cut -f2  | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv

# >>> Reads mapping to tRNA:
grep total \$TEMP/*.tRNA.Summary.txt | cut -f2 -d: | cut -f1 -d\  | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv

# >>> Reads mapping to rRNA:
grep total \$TEMP/*.rRNA.Summary.txt | cut -f2 -d: | cut -f1 -d\  | paste -d, \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv - > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv

# >>> Adding Header:
cat \$FeatC/ReadAssignmentReport_HEAD_\${Mapper}_s\${MapStrand}.csv \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv > tem_\${RandomSeed}; mv tem_\${RandomSeed} \$FeatC/ReadAssignmentReport_\${Mapper}_s\${MapStrand}.csv
rm \$FeatC/ReadAssignmentReport_HEAD_\${Mapper}_s\${MapStrand}.csv

# >>> FeatureCounts Summary: STAR --OR-- tophat2
[[ \$Mapper == "STAR" ]]  && MapResults="Aligned.sortedByCoord.out.bam" || MapResults="/accepted_hits.bam"
sed  's|\t|,|g' \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.txt.summary | \
  sed 's|'\$Map'/||g' |                   \
  sed 's|'\$MapResults'||g' > \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.summary.csv

# >>>>>> Differential Expression Part:
sed -i 's|'\$Map/'||g' \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.txt
sed -i 's|'\$MapResults'||g' \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.txt
head -n2 \$FeatC/\${Mapper}_gene_counts_s\${MapStrand}.txt | tail -n1 | cut -f7- | tr '\t' ',' > \$TEMP/HeadNames.txt

# >>>>>> Run Script for Differential Expression Analysis:
[[ \$DEGanal == 1 ]]   &&  \$R CMD BATCH $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.R $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.Rout

# >>>>>> Run Group Quality control analysis for all and each sample:
[[ \$runQC == "1" ]]  && cd \$QC/geneBody_coverage 
[[ \$runQC == "1" ]]  && /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/geneBody_coverage.py   -i \$Map/MappingNames.txt     -r \$RefSeq -f pdf -o \$QC/geneBody_coverage/Profile  > \$QC/geneBody_coverage/Profile.log           2>&1 # EOF
[[ \$runQC == "1" ]]  && cd \$QC/Transcript_Integrity_Number
[[ \$runQC == "1" ]]  && /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/RSeQC_py3/bin/tin.py                 -i \$Map/MappingNames.txt     -r \$RefSeq                                           > \$QC/Transcript_Integrity_Number/Profile.txt 2>&1 # EOF
[[ \$runQC == "1" ]]  && cat \$QC/Transcript_Integrity_Number/*summary.txt | sort | uniq | tr '\t' ',' | sed 's|Aligned.sortedByCoord.out.bam||' > Summary.csv


EOF
if [ "$run" = "1" ]; then
#  echo "Running Job"
 qsub $Jobs/RNAseq.v${Ver}.${GenVer}.${Mapper}.s${MapStrand}.sh
fi

exit

EOF>
echo 
echo 
echo 

tail -n+3 /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/Xiaocui/02.CD4pos_CsAIono_Stimulation_Across_Time/01.RNAseq/6_6_19_Xiaocui_He_Truseq/04.FeatureCounts/STAR_gene_counts_s0.txt | cut -f1,2 |       cut -f1 -d\;       > $TEMP/Gintervals.01.txt 
tail -n+3 /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/Xiaocui/02.CD4pos_CsAIono_Stimulation_Across_Time/01.RNAseq/6_6_19_Xiaocui_He_Truseq/04.FeatureCounts/STAR_gene_counts_s0.txt | cut -f3   | rev | cut -f1 -d\; | rev > $TEMP/Gintervals.02.txt 
tail -n+3 /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/Xiaocui/02.CD4pos_CsAIono_Stimulation_Across_Time/01.RNAseq/6_6_19_Xiaocui_He_Truseq/04.FeatureCounts/STAR_gene_counts_s0.txt | cut -f4   | rev | cut -f1 -d\; | rev > $TEMP/Gintervals.03.txt 
paste $TEMP/Gintervals.0[123].txt > $TEMP/GenomicIntervals.txt

echo 'sampleid timepoint group
MCsAIono.2h 2h 2
MCsAIono.6h 6h 3
MCsAIono.18h 18h 1
MCsAIono.2h 2h 2
MCsAIono.6h 6h 3
MCsAIono.18h 18h 1
MIono.2h 2h 5
MIono.6h 6h 6
MIono.18h 18h 4
MIono.2h 2h 5
MIono.6h 6h 6
MIono.18h 18h 4
MRest.2h 2h 8
MRest.6h 6h 9
MRest.18h 18h 7
MRest.2h 2h 8
MRest.6h 6h 9
MRest.18h 18h 7' |tr ' ' '\t'  > $TEMP/experiment.txt

tail -n+2 /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/Xiaocui/02.CD4pos_CsAIono_Stimulation_Across_Time/01.RNAseq/6_6_19_Xiaocui_He_Truseq/04.FeatureCounts/STAR_gene_counts_s0.txt | cut -f7- >  $TEMP/counts.txt 

load("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/Xiaocui/02.CD4pos_CsAIono_Stimulation_Across_Time/01.RNAseq/6_6_19_Xiaocui_He_Truseq/Jobs/RNAseq.mm10.STAR.s0.v1.RData")
Cnt = as.matrix(read.table("/mnt/BioScratch/edahi/counts.txt",header=TRUE,stringsAsFactors=FALSE))
Xpr = read.table("/mnt/BioScratch/edahi/experiment.txt", header=TRUE)
Gnt = read.table("/mnt/BioScratch/edahi/GenomicIntervals.txt")
rownames(Cnt) = Gnt[,1]
colnames(Gnt) = c("id","chr","start","end")
tca <- TCA(design = Xpr, genomicFeature = Gnt, counts = Cnt)
tca <- DBanalysis(tca, filter.type = "raw", filter.value = 10, samplePassfilter = 3)
tca <- DBanalysis(tca)

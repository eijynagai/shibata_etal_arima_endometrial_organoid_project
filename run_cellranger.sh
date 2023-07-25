cellranger=/work/CellRanger/cellranger-current/cellranger
reference=/work/Database/10XGenomics/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
wd=/work/nagai/arima_endo_organoid/data/Arima

cd $wd

mkdir -p log

# scRNA-seq
for dir in `ls -d fastq/* | grep -v md5`
do
    label=`basename $dir`
    echo $label
    $cellranger count --id=$label --transcriptome=$reference --fastqs=$dir --sample=$label >log/$label
done

#sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/shortcake.1.3.1.sif"
sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/rnakato_singlecell_jupyter.2022.03.sif"
gtf=/work/Database/10XGenomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf
repeat=/work/Database/10XGenomics/repeat_for_velocyto/hg38_rmsk.gtf
wd=/work3/SingleCell/Arima/data/ari_org

cd $wd

mkdir -p log

#for dir in `ls -d fastq/* | grep -v md5`
for dir in EMO6_hor  Pre1_con  Pre1_hor
do
   # get the directory name	
   #label=`basename $dir`
   #echo $label
   $sing velocyto run10x -m $repeat -@ 24 $dir $gtf #>& log/velocyto.$label.log &
done

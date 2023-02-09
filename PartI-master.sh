#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p mhgcp
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220611_bwt_germ_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220611_bwt_germ_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

BCLDIR="/storage/goodell/bcl/220601_Chunwei_RNA_Trx/Files"
FASTQDIR="/storage/goodell/projects/chunweic/220519_ATACseq_LSK_1"
OUTDIR="/storage/goodell/projects/chunweic/220519_ATACseq_LSK_1/Align"
HOMEDIR="/storage/goodell/home/chunweic"

bcl2fastq -R $BCLDIR -i "$BCLDIR/Data/Intensities/BaseCalls/" -o $FASTQDIR --sample-sheet "$PROJECTDIR/SampleSheet_220601.csv" --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --loading-threads 8 -p 16 --tiles s_[1]

cp -r "$FASTQDIR/220601" $PROJECTDIR
fastqc $PROJECTDIR/*gz
gunzip $PROJECTDIR/*gz
multiqc $PROJECTDIR/*zip -o $PROJECTDIR/

for FILE in $PROJECTDIR/*R1.fastq ;
do
  trimmomatic PE -threads 1 -phred33 \
  $FASTQDIR/Trx_WT1_R1.fastq $FASTQDIR/Trx_WT1_R2.fastq \
  $OUTDIR/Trx_WT1_R1_paired.fastq $OUTDIR/Trx_WT1_R1_unpaired.fastq $OUTDIR/Trx_WT1_R2_paired.fastq $OUTDIR/Trx_WT1_R2_unpaired.fastq \
  ILLUMINACLIP:$HOMEDIR/adapter_trim/NexteraPE-PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

  bowtie2 --dovetail --phred33 -x $HOMEDIR/mm10/BowtieIndex/mm10/mm10 -1 $OUTDIR/Trx_WT1_R1_paired.fastq -2 $OUTDIR/Trx_WT1_R2_paired.fastq -S $OUTDIR/Trx_WT1.sam

  samtools view -S -b $OUTDIR/Trx_WT1.sam > $OUTDIR/Trx_WT1.bam
  samtools sort $OUTDIR/Trx_WT1.bam $OUTDIR/Trx_WT1_sort
  samtools index $OUTDIR/Trx_WT1_sort.bam
done

rm $PROJECTDIR/*fastq $$PROJECTDIR/*zip $$PROJECTDIR/*out.bamn

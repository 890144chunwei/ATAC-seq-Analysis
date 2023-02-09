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

rm $PROJECTDIR/*fastq $$PROJECTDIR/*zip $$PROJECTDIR/*out.bam

#samtools rmdup -S $OUTDIR/Trx_WT3.bam $OUTDIR/Trx_WT3_drm.bam
#samtools sort $OUTDIR/Trx_WT3_drm.bam $OUTDIR/Trx_WT3_drm_sort
#samtools index $OUTDIR/Trx_WT3_drm_sort.bam
#alignmentSieve --numberOfProcessors 8 --ATACshift --bam $OUTDIR/Trx_WT3_drm_sort.bam -o $OUTDIR/Trx_WT3_drm_shift.bam
#samtools sort $OUTDIR/Trx_WT3_drm_shift.bam $OUTDIR/Trx_WT3_drm_shift_sort
#samtools index $OUTDIR/Trx_WT3_drm_shift_sort.bam
#bedtools bamtobed -i $OUTDIR/Trx_WT3_drm_shift_sort.bam > $OUTDIR/Trx_WT3.bed
#bedtools intersect -v -a $OUTDIR/Trx_WT3.bed -b $HOMEDIR/mm10-blacklist.v2.bed > $OUTDIR/Trx_WT3_bl.bed
bamCoverage -b $OUTDIR/Trx_WT3_drm_shift_sort.bam -bl $HOMEDIR/mm10-blacklist.v2.bed -o $OUTDIR/Trx_WT3_cov.bw --normalizeUsing RPKM

#macs3 callpeak -g mm -f BAMPE -t $OUTDIR/Trx_WT1_drm_shift_sort.bam $OUTDIR/Trx_WT2_drm_shift_sort.bam $OUTDIR/Trx_WT3_drm_shift_sort.bam -q 0.01 -n $OUTDIR/Trx_WT --nomodel --keep-dup=all --call-summits
#annotatePeaks.pl $OUTDIR/Trx_WT_summits.bed mm10 -gtf $HOMEDIR/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/Trx_WT.txt

#awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/Trx_WT_peaks.narrowPeak > $OUTDIR/Trx_WT.saf
#featureCounts -T 7 -p -F SAF -a $OUTDIR/Trx_WT.saf -o $OUTDIR/Trx_WT1_fc.txt $OUTDIR/Trx_WT1_drm_shift_sort.bam
#featureCounts -T 7 -p -F SAF -a $OUTDIR/Trx_WT.saf -o $OUTDIR/Trx_WT2_fc.txt $OUTDIR/Trx_WT2_drm_shift_sort.bam
#featureCounts -T 7 -p -F SAF -a $OUTDIR/Trx_WT.saf -o $OUTDIR/Trx_WT3_fc.txt $OUTDIR/Trx_WT3_drm_shift_sort.bam

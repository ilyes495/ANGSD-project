#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH -w buddy.pbtech
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=ProcessAll
#SBATCH --time=48:00:00 # HH/MM/SS
#SBATCH --mem=36G # memory requested, units available: K,M,G,T
outFile=ProcessAll_${SLURM_JOB_ID}.txt
echo "SLURM_SUBMIT_DIR:" ${SLURM_SUBMIT_DIR}  >> ${outFile}
echo "Starting at:" `date` >> ${outFile}
echo "This is job #:" $SLURM_JOB_ID >> ${outFile}
echo "Running on node:" `hostname` >> ${outFile}
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> ${outFile}
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> ${outFile}

spack load -r fastqc
spack load -r trimgalore
spack load star@2.7.0e
spack load samtools@1.9% gcc@6.3.0
spack load subread


#mkdir /scratchLocal/ilb4001
tmp=/scratchLocal/ilb4001
fatqDir=${tmp}/raw_data
#mkdir ${fatqDir}

#echo "Copying fastq files to" ${fatqDir} >> ${outFile}
#cp $SLURM_SUBMIT_DIR/raw_data/SRR67*.fastq.gz ${fatqDir}

#mkdir ${tmp}/fastqc
QCDir=${tmp}/fastqc

#echo "Start fastqc...." >> ${outFile}
#fastqc $(ls ${fatqDir}/SRR67*.fastq.gz) --extract -o ${QCDir} |$ tee ${outFile}

#mkdir ${tmp}/fastq_trimmed
TrimDir=${tmp}/fastq_trimmed
#echo "Start TrimGalore" >> ${outFile}
# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}" --paired  --stringency 5 -o ${TrimDir}  $(ls ${fatqDir}/SRR67*.fastq.gz)
# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}" --paired  --stringency 5 -o ${TrimDir}  $(printf "${tmp}/%s\n" $(cat ${tmp}/raw_fastq.txt))

# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}" --paired  --stringency 5 -o ${TrimDir}  $(ls ${fatqDir}/SRR6703973*.fastq.gz)
# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}"   --stringency 5 -o ${TrimDir}  $(ls ${fatqDir}/SRR6703972*.fastq.gz)
# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}"   --stringency 5 -o ${TrimDir}  $(ls ${fatqDir}/SRR6703979*.fastq.gz)
# trim_galore --gzip --illumina  --fastqc --fastqc_args "--outdir ${QCDir}"   --stringency 5 -o ${TrimDir}  $(ls ${fatqDir}/SRR6703981*.fastq.gz)



#mkdir ${tmp}/genome
GenomeDir=${tmp}/genome

#echo "Copying Genome files" >> ${outFile}
#cp -r $SLURM_SUBMIT_DIR/genome/hg38_STARindex ${GenomeDir}
#cp $SLURM_SUBMIT_DIR/genome/hg38.ensGene.gtf ${GenomeDir}


AlignmDir=${tmp}/alignments
#mkdir ${AlignmDir}
#mkdir ${AlignmDir}/bamqc/
#for sra in $(cat $SLURM_SUBMIT_DIR/raw_data/SRR_Acc_List_v2.txt); do
# for sra in SRR6703973; do
#     sampleName=$(basename ${TrimDir}/${sra}*_1.fq.gz | cut -d '.' -f 1)
    
#     echo "Run STAR on " ${sampleName} >> ${outFile}
#     STAR --runMode alignReads \
#     --runThreadN 16 \
#     --readFilesIn ${TrimDir}/${sra}*_val_1.fq.gz  ${TrimDir}/${sra}*_val_2.fq.gz \
#     --readFilesCommand zcat \
#     --genomeDir ${GenomeDir}/hg38_STARindex \
#     --outSAMtype BAM SortedByCoordinate \
#     --outFileNamePrefix ${AlignmDir}/${sampleName}.
    
#     echo "Run samtools index on " ${sampleName} >> ${outFile}
#     samtools index ${AlignmDir}/${sampleName}.Aligned.sortedByCoord.out.bam
    
#     echo "Run bamqc on " ${sampleName} >> ${outFile}
#     /softlib/apps/EL7/BamQC/bin/bamqc ${AlignmDir}/${sampleName}.Aligned.sortedByCoord.out.bam -o ${AlignmDir}/bamqc/
# done

# for sra in SRR6703972 SRR6703979 SRR6703981; do
#     sampleName=$(basename ${TrimDir}/${sra}*_1.fq.gz | cut -d '.' -f 1)
    
#     echo "Run STAR on " ${sampleName} >> ${outFile}
#     STAR --runMode alignReads \
#     --runThreadN 16 \
#     --readFilesIn ${TrimDir}/${sra}*_val_1.fq.gz \
#     --readFilesCommand zcat \
#     --genomeDir ${GenomeDir}/hg38_STARindex \
#     --outSAMtype BAM SortedByCoordinate \
#     --outFileNamePrefix ${AlignmDir}/${sampleName}.
    
#     echo "Run samtools index on " ${sampleName} >> ${outFile}
#     samtools index ${AlignmDir}/${sampleName}.Aligned.sortedByCoord.out.bam
    
#     echo "Run bamqc on " ${sampleName} >> ${outFile}
#     /softlib/apps/EL7/BamQC/bin/bamqc ${AlignmDir}/${sampleName}.Aligned.sortedByCoord.out.bam -o ${AlignmDir}/bamqc/
# done


GTF=${GenomeDir}/hg38.ensGene.gtf
COUNTFILE=${tmp}/geneCounts.txt 
# Run featureCounts
featureCounts -a ${GTF} \
            -o ${COUNTFILE} \
            -t 'exon' \
            -M \
            --verbose \
            ${AlignmDir}/*.bam



echo "coying the results back...." >> ${outFile}
#cp -r ${AlignmDir}/* ${SLURM_SUBMIT_DIR}/alignments

#mkdir ${SLURM_SUBMIT_DIR}/qc
#cp -r ${QCDir}/* ${SLURM_SUBMIT_DIR}/qc

cp ${tmp}/geneCounts* ${SLURM_SUBMIT_DIR}

#cp -r ${TrimDir} ${SLURM_SUBMIT_DIR}



# echo "remove tmp directory" >> ${outFile}
# rm -r ${tmp}
echo "Done!" >> ${outFile}

exit
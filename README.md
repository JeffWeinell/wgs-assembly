# wgs assembly


###### short reads, iterative reference mapping




pu2b.sh

```
# name of sample
NAME=

# file with settings
SETTINGS_PATH=settings.sh

# load settings
source "$SETTINGS_PATH" "$NAME"


# Iteration 2 Variant calling
# reference genome scaffold names from shortest to longest scaffold length
REFSEQNAMES_02_ALL=$(cat "${FNA_FAI_02}" | awk '{print $2"\t"$1}' | sort -V -r | awk '{print $2}')
NUM_REFSCAFFS_02=$(echo "$REFSEQNAMES_02_ALL" | wc -l)

# Exit if output file from previous step doesn't exist
[ ! -f ${MERGED_RG_MKDUP_02_BAM_PATH} ] && exit 1

# If no merged variant VCFs exist, then do per-scaffold variant calling
if [ ! -f ${FILTER_VCF_02} ]
then
     # iteration 2 loop
     for i in $(seq 1 $NUM_REFSCAFFS_02)
     do
          # name of jth=ith reference scaffold
          REFSEQNAMEj=$(echo "$REFSEQNAMES_02_ALL" | awk -v i=$i 'NR==i{print}')
          # $REFSEQNAMEj scaffold length
          #REFSEQNAMEj_LENGTH=$(cat "${FNA_FAI_02}" | grep $REFSEQNAMEj | awk '{print $2}')
          REFSEQNAMEj_LENGTH=$(awk -v x=$REFSEQNAMEj '{if ($1 == x) print $2}' "${FNA_FAI_02}")
          [ "$REFSEQNAMEj_LENGTH" -lt 10000000 ] && NATIVE_PAIR_HMM_THREADS="$MINPROCS"
          [ "$REFSEQNAMEj_LENGTH" -gt 9999999 ] && [ "$REFSEQNAMEj_LENGTH" -lt 100000000 ] && NATIVE_PAIR_HMM_THREADS="$MIDPROCS"
          [ "$REFSEQNAMEj_LENGTH" -gt 99999999 ] && NATIVE_PAIR_HMM_THREADS="$MAXPROCS"

          GATK_PROCS=$NATIVE_PAIR_HMM_THREADS
          # Rough estimate of number of seconds to run 'gatk HaplotypeCaller'
          GATK_TIME_SEC=$(awk -v REFSEQNAMEj_LENGTH=$REFSEQNAMEj_LENGTH -v GATK_PROCS=$GATK_PROCS 'BEGIN {a=0.0002521; b=REFSEQNAMEj_LENGTH; c=249; d=60; m=GATK_PROCS; print (((a*b)+c)*(d/m)) }' | awk '{printf "%.0f\n", $1}')

          # Require at least one hour per job
          [ "$GATK_TIME_SEC" -lt 3600 ] && GATK_TIME_SEC=3600
          GATK_TIME_SETj="00:00:"$GATK_TIME_SEC
          
          # Where to save pseudo-haplotypes inferred for reads mapped to reference scaffold $REFSEQNAMEj
          GVCFSCAFFj=${OUTDIR}"/iter-02/vcf/gvcf-scaff/"${REFSEQNAMEj}"-iter-02.gvcf.gz"

          # Where to save pseudo-genotypes inferred for reads mapped to reference scaffold $REFSEQNAMEj
          VCFSCAFFj=${OUTDIR}"/iter-02/vcf/gvcf-scaff/"${REFSEQNAMEj}"-iter-02.vcf.gz"

          # Where to save site-filtered (as annotations) pseudo-haplotypes aligned to reference scaffold $REFSEQNAMEj
          VCFSCAFFj_FILTER=${OUTDIR}"/iter-02/vcf/gvcf-scaff/"${REFSEQNAMEj}"-iter-02-filter.vcf.gz"
          
          # Where to save 'gatk HaplotypeCaller' logfile for $REFSEQNAMEj
          HAPCALL_02_LOGj="gatk-haplotypcaller-"${REFSEQNAMEj}"-iter-02.log"

          # Where to save 'gatk GenotypeGVCFs' logfile for $REFSEQNAMEj
          GENOTYPE_02_LOGj="gatk-genotypegvcfs-"${REFSEQNAMEj}"-iter-02.log"

          # Where to save 'bcftools filter' logfile for $REFSEQNAMEj
          FILTERVCF_02_LOGj="bcftools-filter-"${REFSEQNAMEj}"-iter-02.log"

          ### Start 'gatk HaplotypeCaller' job for every $REFSEQNAMEj when 'picard MarkDuplicates' job is complete
          if [ ! -f ${GVCFSCAFFj} ] || [ ! -f ${GVCFSCAFFj}".tbi" ]
          then
               [ -f ${GVCFSCAFFj} ] && rm ${GVCFSCAFFj}
               [ -f ${GVCFSCAFFj}".tbi" ] && rm ${GVCFSCAFFj}".tbi"
               HAPCALL_ID2=$(sbatch --tasks-per-node=$GATK_PROCS --mem=10gb --time=$GATK_TIME_SETj --chdir=${VCFLOGSDIR_02} -o ${HAPCALL_02_LOGj} -e ${HAPCALL_02_LOGj} $HAPLOTYPECALLER_02_SHPATH  ${NAME} ${SETTINGS_PATH} ${NATIVE_PAIR_HMM_THREADS} ${REFSEQNAMEj} ${GVCFSCAFFj} | awk '{print $4}')
               echo -e ${HAPCALL_ID2}"\tgatk_HaplotypeCaller_02 for "${REFSEQNAMEj}"\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 3
          else
               HAPCALL_ID2=$(cat $LOGOFLOGS  | grep "gatk_HaplotypeCaller_02 for "${REFSEQNAMEj} | tail -n 1 | awk '{print $1}')
          fi

          ### Start 'gatk GenotypeGVCFs' job for every $REFSEQNAMEj
          if [ ! -f ${VCFSCAFFj} ] || [ ! -f ${VCFSCAFFj}".tbi" ]
          then
               [ -f ${VCFSCAFFj} ] && rm ${VCFSCAFFj}
               [ -f ${VCFSCAFFj}".tbi" ] && rm ${VCFSCAFFj}".tbi"
               GENOTYPE_ID2=$(sbatch --tasks-per-node=$GATK_PROCS --mem=10gb --time=06:00:00 --chdir=${VCFLOGSDIR_02} -o ${GENOTYPE_02_LOGj} -e ${GENOTYPE_02_LOGj} --dependency=afterok:${HAPCALL_ID2} $GENOTYPEGVCF_02_SHPATH  ${NAME} ${SETTINGS_PATH} ${REFSEQNAMEj} ${GVCFSCAFFj} ${VCFSCAFFj} ${HAPCALL_ID2} | awk '{print $4}')
               echo -e ${GENOTYPE_ID2}"\tgatk_GenotypeGVCFs_02 for "${REFSEQNAMEj}"\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 3
          else
               GENOTYPE_ID2=$(cat $LOGOFLOGS  | grep "gatk_GenotypeGVCFs_02 for "${REFSEQNAMEj} | tail -n 1 | awk '{print $1}')
          fi

          ### start 'bcftools filter' job for $REFSEQNAMEj when 'gatk HaplotypeCaller' job for $REFSEQNAMEj is complete
          if [ ! -f ${VCFSCAFFj_FILTER} ]
          then
               FILTERVCF_ID2=$(sbatch --tasks-per-node=$FILTER_PROCS --mem=500mb --time=01:00:00 --chdir=${VCFLOGSDIR_02} -o ${FILTERVCF_02_LOGj} -e ${FILTERVCF_02_LOGj} --dependency=afterok:${GENOTYPE_ID2} $BCFTOOLSFILTER_02_SHPATH ${NAME} ${SETTINGS_PATH} ${REFSEQNAMEj} ${VCFSCAFFj} ${VCFSCAFFj_FILTER} ${GENOTYPE_ID2} | awk '{print $4}')
               echo -e ${FILTERVCF_ID2}"\tbcftools_filter_02 for "${REFSEQNAMEj}"\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 3
          else
               FILTERVCF_ID2=$(cat $LOGOFLOGS  | grep "bcftools_filter_02 for "${REFSEQNAMEj} | tail -n 1 | awk '{print $1}')
          fi
     
     done

     # List of slurm job IDs for every 'bcftools filter' job (Colon-separated slurm job IDs. Variable that you want to use in the post-loop sbatch dependency)
     BCFTOOLS_FILTER_02_LOGLINES=$(awk -F '\t' '{if ($1 > 1 && $2 ~ "bcftools_filter_02" && $3 == "SUBMITTED" ) print}' "${LOGOFLOGS}")
     FILTERVCF_02_IDS_ALL=$(echo $(echo "$BCFTOOLS_FILTER_02_LOGLINES" | sort -r -k 1,1 | awk -F '\t' '!seen[$2]++' | sort -k 1,1 | awk '{print $1}') | awk 'gsub(" ",":")1')

     ### Create arguments file for 'gatk GatherVCFs'
     if [ ! -f ${ARGSFILE_02} ]
     then
          CREATEARGS_02_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=00:15:00 --chdir=$OUTDIR"/iter-02/logs/" -o "create_argfile_02.log" -e "create_argfile_02.log" --dependency=afterok:${FILTERVCF_02_IDS_ALL} $CREATEARGS_02_SHPATH ${NAME} ${SETTINGS_PATH} ${FILTERVCF_02_IDS_ALL} | awk '{print $4}')
          echo -e ${CREATEARGS_02_ID}"\tcreate_argfile_02\tSUBMITTED\t-----" >> ${LOGOFLOGS}
          sleep 1
     else
          # Get job id from log of logs
          CREATEARGS_02_ID=$(cat $LOGOFLOGS  | grep "create_argfile_02" | tail -n 1 | awk '{print $1}')
     fi

     GATHERVCFS_02_ID=$(sbatch --tasks-per-node=1 --mem=12gb --time=01:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "GatherVCFs_02.log" -e "GatherVCFs_02.log" --dependency=afterok:${CREATEARGS_02_ID} $GATHERVCFS_02_SHPATH ${NAME} ${SETTINGS_PATH} ${CREATEARGS_02_ID} | awk '{print $4}')
     echo -e ${GATHERVCFS_02_ID}"\tGatherVCFs_02\tSUBMITTED\t-----" >> ${LOGOFLOGS}
     sleep 1

else
     # Get/set variable with completed job ids from log of logs
     GATHERVCFS_02_ID=$(cat $LOGOFLOGS  | grep "GatherVCFs_02" | tail -n 1 | awk '{print $1}')
fi

```

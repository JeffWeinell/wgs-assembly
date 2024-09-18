# wgs assembly


#### short reads, iterative reference mapping

settings.sh
```

```


pu1a.sh

```
NAME=
SETTINGS_PATH=

# Load environment and assembly settings
source ${SETTINGS_PATH} ${NAME}

#################################################
############### Start iteration 1 ###############
#################################################

# check/create copy of reference fasta and create dictionary/index files
if [ ! -d "$REFDIR_01" ]
then
     CPREF_01_ID=$(sbatch --tasks-per-node=10 --mem=10gb --time=04:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "tcp_ref_01.log" -e "tcp_ref_01.log" $CPREF_01_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
     echo -e ${CPREF_01_ID}"\tcp_ref_01\tSUBMITTED\t---" >> $LOGOFLOGS
     sleep 1
else
     # Get job id from log of logs
     if grep -q "cp_ref_01" $LOGOFLOGS
     then
          CPREF_01_ID=$(cat $LOGOFLOGS | grep "cp_ref_01" | tail -n 1 | awk '{print $1}')
     else
          CPREF_01_ID=${SLURM_JOB_ID} && echo -e ${CPREF_01_ID}"\tcp_ref_01\tSUBMITTED\t---" >> $LOGOFLOGS
     fi
fi

# check for bamfile with reads mapped to reference and duplicate reads marked
if [ ! -f $OUTDIR"/iter-01/bam/merged-rg-mkdup-iter-01.bam.gz" ]
then
     # check for bamfile with SE and PE reads mapped to reference
     if [ ! -f $OUTDIR"/iter-01/bam/merged-iter-01.bam.gz" ]
     then
          # SE reads mapping
          if [ ! -f $OUTDIR"/iter-01/bam/se-iter-01.bam.gz" ]
          then
               # 'bwa mem ' for SE reads mapping to reference genome scaffolds
               BWA_MEM_SE_ID=$(sbatch --tasks-per-node=${BWA_THREADS} --mem=100gb --time=200:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "bwa_mem_se_01.log"  -e "bwa_mem_se_01.log" --dependency=afterok:${CPREF_01_ID} $BWA_MEM_SE_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
               echo -e ${BWA_MEM_SE_ID}"\tbwa_mem_se_01\tSUBMITTED\t---" >> $LOGOFLOGS
               sleep 1
          else 
               # Get job id from log of logs
               BWA_MEM_SE_ID=$(cat $LOGOFLOGS | grep "bwa_mem_se_01" | tail -n 1 | awk '{print $1}')
          fi

          # PE reads mapping
          if [ ! -f $OUTDIR"/iter-01/bam/pe-iter-01.bam.gz" ]
          then
               # 'bwa mem ' for PE reads mapping to reference genome scaffolds
               BWA_MEM_PE_ID=$(sbatch --tasks-per-node=${BWA_THREADS} --mem=100gb --time=200:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "bwa_mem_pe_01.log"  -e "bwa_mem_pe_01.log" --dependency=afterok:${CPREF_01_ID} $BWA_MEM_PE_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
               echo -e ${BWA_MEM_PE_ID}"\tbwa_mem_pe_01\tSUBMITTED\t---" >> $LOGOFLOGS
               sleep 1
          else
               # Get job id from log of logs
               BWA_MEM_PE_ID=$(cat $LOGOFLOGS | grep "bwa_mem_pe_01" | tail -n 1 | awk '{print $1}')
          fi
          
          # Merge SE and PE bams
          MERGESAM_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=05:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "picard_mergesam_01.log"  -e "picard_mergesam_01.log" --dependency=afterok:${BWA_MEM_SE_ID}:${BWA_MEM_PE_ID} $MERGESAM_01_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
          echo -e ${MERGESAM_01_ID}"\tpicard_mergesam_01\tSUBMITTED\t---" >> $LOGOFLOGS
          sleep 1
     else
          # Get job id from log of logs
          MERGESAM_01_ID=$(cat $LOGOFLOGS | grep "picard_mergesam_01" | tail -n 1 | awk '{print $1}')
          # MERGESAM_01_ID="$CPREF_01_ID"
     fi
     
     # Mark duplicate reads
     MARKDUP_ID=$(sbatch --tasks-per-node=10 --mem=100gb --time=5:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "mark_duplicates_01.log"  -e "mark_duplicates_01.log" --dependency=afterok:${MERGESAM_01_ID} $MARKDUP_01_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
     echo -e ${MARKDUP_ID}"\tmark_duplicates_01\tSUBMITTED\t---" >> $LOGOFLOGS
     sleep 1
     
else
     # Get job id from log of logs
     MARKDUP_ID=$(cat $LOGOFLOGS | grep "mark_duplicates_01" | tail -n 1 | awk '{print $1}')
     # MARKDUP_ID="$CPREF_01_ID"
fi

######################################
# sbatch ------> i1 haplotyping swarm submitter
######################################

PU1B_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=6:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "PU1B.log"  -e "PU1B.log" --dependency=afterok:${MARKDUP_ID} $PU1B_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
echo -e ${PU1B_ID}"\tpu1b\tSUBMITTED\t---" >> $LOGOFLOGS

exit 0
```

pu1b.sh

```
NAME=
SETTINGS_PATH=


# exit with error if NAME argument empty
[ ! ${NAME} ] && echo "NAME argument empty" && exit 1

# Load environment and assembly settings
source ${SETTINGS_PATH} ${NAME}

# check that all previous jobs completed successfully
if grep -q "mark_duplicates_01" $LOGOFLOGS
then
     MARKDUP_ID=$(cat $LOGOFLOGS | grep "mark_duplicates_01" | tail -n 1 | awk '{print $1}')
else
     exit 1
fi

# Reference genome scaffold names from shortest to longest scaffold length
REFSEQNAMES_ALL=$(cat "${FNA_FAI}" | awk '{print $2"\t"$1}' | sort -V -r | awk '{print $2}')
NUM_REFSCAFFS=$(echo "$REFSEQNAMES_ALL" | wc -l)
# If no merged variant VCFs exist, then check/do per-scaffold variant calling
if [ ! -f ${FILTER_INTERMEDIATE_VCF_01} ] && [ ! -f ${FILTER_INTERMEDIATE_SNPS_VCF_01} ]
then
     # iteration 1 loop
     for i in $(seq 1 $NUM_REFSCAFFS)
     do
          # name of jth=ith reference scaffold
          REFSEQNAMEj=$(echo "$REFSEQNAMES_ALL" | awk -v i=$i 'NR==i{print}')
          # $REFSEQNAMEj scaffold length
          #REFSEQNAMEj_LENGTH=$(cat "${FNA_FAI}" | grep $REFSEQNAMEj | awk '{print $2}')
          REFSEQNAMEj_LENGTH=$(awk -v x=$REFSEQNAMEj '{if ($1 == x) print $2}' "${FNA_FAI}")
          # How many threads should a native pairHMM implementation use (gatk default = 4)
          [ "$REFSEQNAMEj_LENGTH" -lt 10000000 ] && NATIVE_PAIR_HMM_THREADS="$MINPROCS"
          [ "$REFSEQNAMEj_LENGTH" -gt 9999999 ] && [ "$REFSEQNAMEj_LENGTH" -lt 100000000 ] && NATIVE_PAIR_HMM_THREADS="$MIDPROCS"
          [ "$REFSEQNAMEj_LENGTH" -gt 99999999 ] && NATIVE_PAIR_HMM_THREADS="$MAXPROCS"
          GATK_PROCS=$NATIVE_PAIR_HMM_THREADS
          # Number of seconds to run 'gatk HaplotypeCaller' on the cluster (based on linear regression of scaffold length versus time to complete). Other factors will matter too of course. Increase value of 'd' on the next line any jobs run out of time.
          GATK_TIME_SEC=$(awk -v REFSEQNAMEj_LENGTH=$REFSEQNAMEj_LENGTH -v GATK_PROCS=$GATK_PROCS 'BEGIN {a=0.0002521; b=REFSEQNAMEj_LENGTH; c=249; d=60; m=GATK_PROCS; print (((a*b)+c)*(d/m)) }' | awk '{printf "%.0f\n", $1}')
          # allocate at least one hour for every job
          [ "$GATK_TIME_SEC" -lt 3600 ] && GATK_TIME_SEC=3600
          GATK_TIME_SETj="00:00:"$GATK_TIME_SEC
          # Where to save pseudo-haplotypes inferred for reads mapped to reference scaffold $REFSEQNAMEj
          VCFSCAFFj=${OUTDIR}"/iter-01/vcf/vcf-scaff/"${REFSEQNAMEj}"-iter-01.vcf.gz"
          # Name of HaplotypeCaller_01 logfile for scaffold REFSEQNAMEj
          HAPCALL_LOGj="gatk-haplotypcaller-"${REFSEQNAMEj}"-iter-01.log"
          FILTERVCF_LOGj="bcftools-filter-"${REFSEQNAMEj}"-iter-01.log"
          # check if scaffold vcf already exists (may switch this to look for a log file instead)
          if [ ! -f ${VCFSCAFFj} ] || [ ! -f ${VCFSCAFFj}".tbi" ]
          then
               [ -f ${VCFSCAFFj} ] && rm ${VCFSCAFFj}
               [ -f ${VCFSCAFFj}".tbi" ] && rm ${VCFSCAFFj}".tbi"
               ### 'gatk HaplotypeCaller' job for scaffold $REFSEQNAMEj
               HAPCALL_ID=$(sbatch --tasks-per-node=$GATK_PROCS --mem=2gb --time=$GATK_TIME_SETj --chdir="$VCFLOGSDIR_01" -o "$HAPCALL_LOGj" -e "$HAPCALL_LOGj" --dependency=afterok:${MARKDUP_ID} $HAPLOTYPECALLER_SHPATH ${NAME} ${SETTINGS_PATH} ${REFSEQNAMEj} ${NATIVE_PAIR_HMM_THREADS} ${VCFSCAFFj} | awk '{print $4}')
               echo -e ${HAPCALL_ID}"\tgatk_HaplotypeCaller_01 for "${REFSEQNAMEj}"\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 3
          else
               # Get job id from log of logs
               HAPCALL_ID=$(cat $LOGOFLOGS  | grep "gatk_HaplotypeCaller_01 for "${REFSEQNAMEj} | tail -n 1 | awk '{print $1}')
          fi
          # Where to save site-filtered pseudo-haplotypes aligned to reference scaffold $REFSEQNAMEj
          VCFSCAFFj_FILTER=${OUTDIR}"/iter-01/vcf/vcf-scaff/"${REFSEQNAMEj}"-iter-01-filter.vcf.gz"
          if [ ! -f ${VCFSCAFFj_FILTER} ]
          then
               ### start 'bcftools filter' job for $REFSEQNAMEj when 'gatk HaplotypeCaller' job for $REFSEQNAMEj is complete
               FILTERVCF_ID=$(sbatch --tasks-per-node=$FILTER_PROCS --mem=500mb --time=02:00:00 --chdir="$VCFLOGSDIR_01" -o "$FILTERVCF_LOGj" -e "$FILTERVCF_LOGj" --dependency=afterok:${HAPCALL_ID} $BCFTOOLSFILTER_SHPATH  ${NAME} ${SETTINGS_PATH} ${REFSEQNAMEj} ${VCFSCAFFj} ${VCFSCAFFj_FILTER} ${HAPCALL_ID} | awk '{print $4}')
               echo -e ${FILTERVCF_ID}"\tbcftools_filter_01 for "${REFSEQNAMEj}"\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 3
          else
               # Get job id from log of logs
               FILTERVCF_ID=$(cat $LOGOFLOGS  | grep "bcftools_filter_01 for "${REFSEQNAMEj} | tail -n 1 | awk '{print $1}')
          fi
     done
     
     # List of slurm job IDs for every 'bcftools filter' job
     #### Colon-separated slurm job IDs. Variable that you want to use in the post-loop sbatch dependency.
     #FILTERVCF_01_IDS_ALL=$(echo $(cat "$LOGOFLOGS" | grep bcftools_filter_01 | grep SUBMITTED | awk '{print $1}') | awk 'gsub(" ",":")1')
     #FILTERVCF_01_IDS_ALL=$(echo $(awk -F '\t' '{if ($1 > 1 && $2 ~ "bcftools_filter_01" && $3 == "SUBMITTED" ) print $1}' "${LOGOFLOGS}") | awk 'gsub(" ",":")1')

     BCFTOOLS_FILTER_01_LOGLINES=$(awk -F '\t' '{if ($1 > 1 && $2 ~ "bcftools_filter_01" && $3 == "SUBMITTED" ) print}' "${LOGOFLOGS}")
     FILTERVCF_01_IDS_ALL=$(echo $(echo "$BCFTOOLS_FILTER_01_LOGLINES" | sort -r -k 1,1 | awk -F '\t' '!seen[$2]++' | sort -k 1,1 | awk '{print $1}') | awk 'gsub(" ",":")1')

     ### Create arguments file for 'gatk GatherVCFs'
     if [ ! -f ${ARGSFILE_01} ]
     then
          CREATEARGS_01_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=00:15:00 --chdir=$OUTDIR"/iter-01/logs/" -o "create_argfile_01.log"  -e "create_argfile_01.log" --dependency=afterok:${FILTERVCF_01_IDS_ALL} $CREATEARGS_01_SHPATH ${NAME} ${SETTINGS_PATH} ${FILTERVCF_01_IDS_ALL} | awk '{print $4}')
          echo -e ${CREATEARGS_01_ID}"\tcreate_argfile_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
          sleep 1
     else
          # Get job id from log of logs
          CREATEARGS_01_ID=$(cat $LOGOFLOGS  | grep "create_argfile_01" | tail -n 1 | awk '{print $1}')
     fi
     
     # Gather iteration 1 filtered pseudoscaffold VCFs
     if [ ! -f ${FILTER_INTERMEDIATE_VCF_01} ]
     then
          GATHERVCFS_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=01:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "GatherVCFs_01.log"  -e "GatherVCFs_01.log" --dependency=afterok:${CREATEARGS_01_ID} $GATHERVCFS_01_SHPATH ${NAME} ${SETTINGS_PATH} ${CREATEARGS_01_ID} | awk '{print $4}')
          echo -e ${GATHERVCFS_01_ID}"\tGatherVCFs_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
          sleep 1
     else
          # Get job id from log of logs
          GATHERVCFS_01_ID=$(cat $LOGOFLOGS  | grep "GatherVCFs_01" | tail -n 1 | awk '{print $1}')
     fi
else
     # Get/set variable with completed job ids from log of logs
     GATHERVCFS_01_ID=$CPREF_01_ID
fi

# Index filtered VCF
if [ ! -f ${FILTER_INTERMEDIATE_VCF_01}".tbi" ] && [ ! -f ${FILTER_INTERMEDIATE_SNPS_VCF_01} ]
then
     TABIX_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=01:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "tabix_01.log"  -e "tabix_01.log" --dependency=afterok:${GATHERVCFS_01_ID} $TABIX_01_SHPATH ${NAME} ${SETTINGS_PATH} ${GATHERVCFS_01_ID} | awk '{print $4}')
     echo -e ${TABIX_01_ID}"\ttabix_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
     sleep 1
else
     # Get job id from log of logs
     TABIX_01_ID=$(cat $LOGOFLOGS  | grep "tabix_01" | tail -n 1 | awk '{print $1}')
     # TABIX_01_ID=$CPREF_01_ID
fi

# Select variants
if [ ! -f ${FILTER_INTERMEDIATE_SNPS_VCF_01} ]
then
     SELECTVARIANTS_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=01:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "gatk_SelectVariants_01.log"  -e "gatk_SelectVariants_01.log" --dependency=afterok:${TABIX_01_ID} $SELECTVARIANTS_01_SHPATH ${NAME} ${SETTINGS_PATH} ${TABIX_01_ID} | awk '{print $4}')
     echo -e ${SELECTVARIANTS_01_ID}"\tgatk_SelectVariants_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
     sleep 1
else
     # Get job id from log of logs
     SELECTVARIANTS_01_ID=$(cat $LOGOFLOGS  | grep "gatk_SelectVariants_01" | tail -n 1 | awk '{print $1}')
     # SELECTVARIANTS_01_ID=$CPREF_01_ID
fi

# Index filtered SNPs VCF
if [ ! -f ${FILTER_INTERMEDIATE_SNPS_VCF_01}".tbi" ]
then
     TABIX_SNPS_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=01:00:00 --chdir=$OUTDIR"/iter-01/logs/" -o "tabix_snps_01.log" -e "tabix_snps_01.log" --dependency=afterok:${SELECTVARIANTS_01_ID} $TABIX_SNPS_01_SHPATH ${NAME} ${SETTINGS_PATH} ${SELECTVARIANTS_01_ID} | awk '{print $4}')
     echo -e ${TABIX_SNPS_01_ID}"\ttabix_snps_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
     sleep 1
else
     # Get job id from log of logs
     TABIX_SNPS_01_ID=$(cat $LOGOFLOGS  | grep "tabix_snps_01" | tail -n 1 | awk '{print $1}')
     # TABIX_SNPS_01_ID=$CPREF_01_ID
fi

# Determine case of first base and if lowercase change to uppercase
# See here for some info for 'bcftools consensus' https://github.com/samtools/bcftools/issues/1150
#     From pseudoit getConsCase() documentation: "bcftools consensus uses the case of the first base in the reference file for injecting variants. This is a hack to make sure it is always upper case."
#     If the first base of the ref were allowed to be lower case, then the entire consensus would be lowercase and subsequent softmasking would be lost (it would actually appear that everything is soft-masked).
#### skipping this because cp_ref_01 copied reference genome to all caps
## if grep -q "first_base_upper_01" $LOGOFLOGS
## then
##      FIRSTBASEUPPER_01_ID=$(cat $LOGOFLOGS  | grep "first_base_upper_01" | tail -n 1 | awk '{print $1}')
## else
##      FIRSTBASEUPPER_01_ID=$(sbatch --tasks-per-node=1 --mem=10gb --time=00:10:00 --chdir=$OUTDIR"/iter-01/logs/" -o "first_base_upper_01.log" -e "first_base_upper_01.log" --dependency=afterok:${TABIX_SNPS_01_ID} $FIRSTBASEUPPER_01_SHPATH ${NAME} ${SETTINGS_PATH} ${TABIX_SNPS_01_ID} | awk '{print $4}')
##      echo -e ${FIRSTBASEUPPER_01_ID}"\tfirst_base_upper_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
##      sleep 1
## fi

# Create iteration 1 consensus fasta, which will be used as the reference for iteration 2
if [ ! -f ${REFPATH_02} ]
then
     # BCFTOOLSCONSENSUS_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=00:30:00 --chdir=$OUTDIR"/iter-01/logs/" -o "bcftools_consensus_01.log" -e "bcftools_consensus_01.log" --dependency=afterok:${FIRSTBASEUPPER_01_ID} $BCFTOOLSCONSENSUS_01_SHPATH ${REFPATH_IN} ${REFPATH_02} ${FILTER_INTERMEDIATE_SNPS_VCF_01} ${FIRSTBASEUPPER_01_ID} ${LOGOFLOGS} | awk '{print $4}')
     BCFTOOLSCONSENSUS_01_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=00:30:00 --chdir=$OUTDIR"/iter-01/logs/" -o "bcftools_consensus_01.log" -e "bcftools_consensus_01.log" --dependency=afterok:${TABIX_SNPS_01_ID} $BCFTOOLSCONSENSUS_01_SHPATH ${NAME} ${SETTINGS_PATH} ${TABIX_SNPS_01_ID} | awk '{print $4}')
     echo -e ${BCFTOOLSCONSENSUS_01_ID}"\tbcftools_consensus_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
     sleep 1
else
     BCFTOOLSCONSENSUS_01_ID=$(cat $LOGOFLOGS  | grep "bcftools_consensus_01" | tail -n 1 | awk '{print $1}')
     # BCFTOOLSCONSENSUS_01_ID=$FIRSTBASEUPPER_01_ID
fi

## # Iteration 1 cleanup (for now nothing is deleted)
## if grep -q "cleanup_01" $LOGOFLOGS
## then
##      CLEANUP_01_ID=$(cat $LOGOFLOGS  | grep "cleanup_01" | tail -n 1 | awk '{print $1}')
## else
##      CLEANUP_01_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=00:05:00 --chdir=$OUTDIR"/iter-01/logs/" -o "cleanup_01.log" -e "cleanup_01.log" --dependency=afterok:${BCFTOOLSCONSENSUS_01_ID} $CLEANUP_01_SHPATH ${OUTDIR} ${KEEPALL} ${BCFTOOLSCONSENSUS_01_ID} ${LOGOFLOGS} | awk '{print $4}')
##      echo -e ${CLEANUP_01_ID}"\tcleanup_01\tSUBMITTED\t-----" >> ${LOGOFLOGS}
##      sleep 1
## fi

############### End of iteration 1 ###############

######################################
# sbatch ------> i2 mapping
######################################

PU2A_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=1:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "PU2A.log"  -e "PU2A.log" --dependency=afterok:${BCFTOOLSCONSENSUS_01_ID} $PU2A_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
echo -e ${PU2A_ID}"\tpu2a\tSUBMITTED\t---" >> $LOGOFLOGS

exit 0

```

pu2a.sh

```
NAME=$1
SETTINGS_PATH=$2

# exit with error if NAME argument empty
[ ! ${NAME} ] && echo "NAME argument empty" && exit 1
## echo $NAME >> $OUTLOG # <--- uncomment after uploading to cluster

# Load environment and assembly settings
source ${SETTINGS_PATH} ${NAME}


#################################################
############### Start iteration 2 ###############
#################################################

# check/create copy of reference fasta and create dictionary/index files
if [ ! -d "$REFDIR_02" ]
then
     CPREF_02_ID=$(sbatch --tasks-per-node=10 --mem=10gb --time=04:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "tcp_ref_02.log" -e "tcp_ref_02.log" $CPREF_02_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
     echo -e ${CPREF_02_ID}"\tcp_ref_02\tSUBMITTED\t---" >> $LOGOFLOGS
     sleep 1
else
     # Get job id from log of logs
     if grep -q "cp_ref_02" $LOGOFLOGS
     then
          CPREF_02_ID=$(cat $LOGOFLOGS | grep "cp_ref_02" | tail -n 1 | awk '{print $1}')
     else
          CPREF_02_ID=${SLURM_JOB_ID} && echo -e ${CPREF_02_ID}"\tcp_ref_02\tSUBMITTED\t---" >> $LOGOFLOGS
     fi
fi

# check for bamfile with reads mapped to reference and duplicate reads marked
if [ ! -f $OUTDIR"/iter-02/bam/merged-rg-mkdup-iter-02.bam.gz" ]
then
     # check for bamfile with SE and PE reads mapped to reference
     if [ ! -f $OUTDIR"/iter-02/bam/merged-iter-02.bam.gz" ]
     then
          # SE reads mapping
          if [ ! -f $OUTDIR"/iter-02/bam/se-iter-02.bam.gz" ]
          then
               BWA_MEM_SE_02_ID=$(sbatch --tasks-per-node=${BWA_THREADS} --mem=100gb --time=200:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "bwa_mem_se_02.log" -e "bwa_mem_se_02.log" --dependency=afterok:${CPREF_02_ID} $BWA_MEM_SE_02_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
               echo -e ${BWA_MEM_SE_02_ID}"\tbwa_mem_se_02\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 1
          else
               # Get job id from log of logs
               BWA_MEM_SE_02_ID=$(cat $LOGOFLOGS | grep "bwa_mem_se_02" | tail -n 1 | awk '{print $1}')
          fi
          
          # PE reads mapping
          if [ ! -f $OUTDIR"/iter-02/bam/pe-iter-02.bam.gz" ]
          then
               BWA_MEM_PE_02_ID=$(sbatch --tasks-per-node=${BWA_THREADS} --mem=100gb --time=200:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "bwa_mem_pe_02.log" -e "bwa_mem_pe_02.log" --dependency=afterok:${CPREF_02_ID} $BWA_MEM_PE_02_SHPATH ${NAME} ${SETTINGS_PATH} ${CREATEDICT_02_ID} ${BWAINDEX_02_ID} ${FAIDX_02_ID} | awk '{print $4}')
               echo -e ${BWA_MEM_PE_02_ID}"\tbwa_mem_pe_02\tSUBMITTED\t-----" >> ${LOGOFLOGS}
               sleep 1
          else
               # Get job id from log of logs
               BWA_MEM_PE_02_ID=$(cat $LOGOFLOGS | grep "bwa_mem_pe_02" | tail -n 1 | awk '{print $1}')
          fi
          
          # Merge SE and PE bams
          MERGESAM_02_ID=$(sbatch --tasks-per-node=1 --mem=100gb --time=10:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "picard_mergesam_02.log" -e "picard_mergesam_02.log" --dependency=afterok:${BWA_MEM_SE_02_ID}:${BWA_MEM_PE_02_ID} $MERGESAM_02_SHPATH ${NAME} ${SETTINGS_PATH} ${BWA_MEM_SE_02_ID} ${BWA_MEM_PE_02_ID} | awk '{print $4}')
          echo -e ${MERGESAM_02_ID}"\tpicard_mergesam_02\tSUBMITTED\t---" >> $LOGOFLOGS
          sleep 1
     else
          MERGESAM_02_ID=$(cat $LOGOFLOGS | grep "picard_mergesam_02" | tail -n 1 | awk '{print $1}')
          ## MERGESAM_02_ID="$CLEANUP_01_ID"
     fi

     # Mark duplicates
     MARKDUP_02_ID=$(sbatch --tasks-per-node=2 --mem=100gb --time=23:00:00 --chdir=$OUTDIR"/iter-02/logs/" -o "mark_duplicates_02.log" -e "mark_duplicates_02.log" --dependency=afterok:${MERGESAM_02_ID} $MARKDUP_02_SHPATH ${NAME} ${SETTINGS_PATH} ${MERGESAM_02_ID} | awk '{print $4}')
     echo -e ${MARKDUP_02_ID}"\tmark_duplicates_02\tSUBMITTED\t---" >> $LOGOFLOGS
     sleep 1
else
     # Get job id from log of logs
     MARKDUP_02_ID=$(cat $LOGOFLOGS | grep "mark_duplicates_02" | tail -n 1 | awk '{print $1}')
     ## MARKDUP_02_ID="$CLEANUP_01_ID"
fi

######################################
# sbatch ------> i2 haplotyping swarm submitter
######################################

# schedule pu2b.sh to start running after other jobs finish successfully
PU2B_ID=$(sbatch --tasks-per-node=1 --mem=100mb --time=6:00:00 --chdir=$OUTDIR"/iter-01/logs/"  -o "PU2B.log"  -e "PU2B.log" --dependency=afterok:${MARKDUP_02_ID} $PU2B_SHPATH ${NAME} ${SETTINGS_PATH} | awk '{print $4}')
echo -e ${PU2B_ID}"\tpu2b\tSUBMITTED\t---" >> $LOGOFLOGS

exit 0
```



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

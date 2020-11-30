#!/bin/env bash

#$ -cwd
#$ -V
#$ -pe smp 4
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail

usage="Usage: $0 -i [input_dir] -o [output_dir]"

while getopts "i:o:" opt; do
  case "${opt}" in
    i )
        in_dir=${OPTARG}
      ;;
    o )
        out_dir=${OPTARG}
      ;;
    \? )
        echo ${usage}
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$in_dir" ]] || [[ -z "$out_dir" ]]; then
    echo $usage
    exit 1
fi

if [[ -z "$SGE_TASK_ID" ]] || [[ "$SGE_TASK_ID" == 'undefined' ]]; then
    echo "Error: Script must be submitted as an array job (qsub -t 1-n)"
    exit 1
fi

mkdir -p ${out_dir} 

readarray -t samples < <(ls -1 ${in_dir}/*bam)
bam=$(basename ${samples[$SGE_TASK_ID-1]#*"/"})
sample=${bam%%"."*}

echo
echo "hostname=$HOSTNAME"
echo "bam file=$bam"
echo "sample=$sample"
echo "in_dir=$in_dir"
echo "out_dir=$out_dir"
echo

cp -rv ${in_dir}/${bam}* ${TMPDIR}/
samtools fastq -f 12 -@ 4 -N -1 ${TMPDIR}/${sample}.filtered_1.fq.gz -2 ${TMPDIR}/${sample}.filtered_2.fq.gz -0 /dev/null ${TMPDIR}/${bam}
cp ${TMPDIR}/*fq.gz ${out_dir}


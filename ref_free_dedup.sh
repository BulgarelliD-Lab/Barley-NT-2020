#!/bin/env bash

#$ -cwd
#$ -V
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -f [fastq_dir] -o [output_dir]"
  echo
  echo "Carries out reference-free deduplification of reads using bbtools clumpify.sh"
  echo
  echo "  -f: Fastq directory. Files should be paired and suffixed with _[12].fq.gz"
  echo "  -o: output directory"
  echo
  exit 1
}

while getopts "f:o:" opt; do
  case "${opt}" in
    f ) 
        fastq_dir=${OPTARG}
    ;;
    o )
        out_dir=${OPTARG}
    ;;
    \? )
        display_help
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$fastq_dir" ]] || [[ -z "$out_dir" ]]; then
    display_help
    exit 1
fi

if [[ -z "$in_dir" ]] || [[ -z "$out_dir" ]]; then
    echo $usage
    exit 1
fi

if [[ -z "$SGE_TASK_ID" ]] || [[ "$SGE_TASK_ID" == 'undefined' ]]; then
	echo "Job must be submitted as an array job"
	exit 1
fi

readarray -t samples < <(ls -1 ${fastq_dir}|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]%%'.'*}

echo
echo "hostname=$HOSTNAME"
echo "sample=$sample"
echo "fastq_dir=${fastq_dir}"
echo "out_dir=${out_dir}"
echo

mkdir -p ${out_dir}

cp -v ${fastq_dir}/${sample}* $TMPDIR

clumpify.sh in=$TMPDIR/${sample}.filtered_1.fq.gz in2=$TMPDIR/${sample}.filtered_2.fq.gz \
    out=$TMPDIR/${sample}.nodup_1.fq.gz out2=$TMPDIR/${sample}.nodup_2.fq.gz dedupe=t

cp -v $TMPDIR/${sample}.nodup* ${out_dir}


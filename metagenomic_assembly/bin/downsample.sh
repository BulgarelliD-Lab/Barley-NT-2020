#!/bin/env bash

#$ -cwd
#$ -V
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -i [input_dir] -o [output_dir]"
  echo
  echo "Carries out reference-free deduplification of reads using bbtools clumpify.sh"
  echo
  echo "  -f: input fastq directory. Files should be paired and suffixed with _[12].fq.gz"
  echo "  -o: output directory"
  echo "  -c: count of reads required per library"
  echo
  exit 1
}

while getopts "i:o:c:" opt; do
  case "${opt}" in
    i )
        in_dir=${OPTARG}
      ;;
    o )
        out_dir=${OPTARG}
      ;;
    c )
        read_count=${OPTARG}
      ;;
    \? )
        echo ${usage}
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$in_dir" ]] || [[ -z "$out_dir" ]] || [[ -z "$read_count" ]]; then
    echo $usage
    exit 1
fi

if [[ -z "$SGE_TASK_ID" ]] || [[ "$SGE_TASK_ID" == 'undefined' ]]; then
    echo "Error: Script must be submitted as an array job (qsub -t 1-n)"
    exit 1
fi

readarray -t samples < <(ls -1 ${in_dir}|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]}

echo
echo "hostname=$HOSTNAME"
echo "Sample=${sample}"
echo "fastq_dir=${fastq_dir}"
echo "out_dir=${out_dir}"
echo "count=${read_count}"
echo

mkdir -p ${out_dir}

cp -v ${in_dir}/${sample}* $TMPDIR
r1=$(ls $TMPDIR/${sample}*1.fq.gz)
r2=$(ls $TMPDIR/${sample}*2.fq.gz)

echo "r1=$r1"
echo "r2=$r2"
echo
reformat.sh in1=${r1} in2=${r2} out1=$TMPDIR/${sample}_downsampled_1.fq.gz \
    out2=$TMPDIR/${sample}_downsampled_2.fq.gz samplereadstarget=${read_count}
cp $TMPDIR/${sample}_downsampled_?.fq.gz ${out_dir}


#!/bin/env bash

#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 16
#$ -mods l_hard h_vmem 48G
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -i [input_dir] -o [output_dir]"
  echo
  echo "Carries out reference-free deduplification of reads using bbtools clumpify.sh"
  echo
  echo "  -f: input fastq directory. Files should be paired and suffixed with _[12].fq.gz"
  echo "  -o: output directory"
  echo
  exit 1
}

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

readarray -t samples < <(ls -1 ${in_dir}|cut -f1 -d_|sort -u)
sample=${samples[$SGE_TASK_ID-1]%%'.'*}

echo "hostname=$HOSTNAME"
echo "sample=$sample"
echo "in_dir=$in_dir"
echo "out_dir=$out_dir"
echo

cp -v ${in_dir}/${sample}*.fq.gz $TMPDIR
r1=$(ls $TMPDIR/${sample}*1.fq.gz)
r2=$(ls $TMPDIR/${sample}*2.fq.gz)

echo "r1=$r1"
echo "r2=$r2"
echo

metaspades.py -t 16 -o ${TMPDIR}/${sample} -1 ${r1} -2 ${r2}
rm $TMPDIR/*fq.gz
cp -rv ${TMPDIR}/${sample} ${out_dir}


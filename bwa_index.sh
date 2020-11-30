#!/bin/env bash

#$ -cwd
#$ -V
#$ -jc long
#$ -mods l_hard mfree 64G
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -r [reference_fasta]"
  echo
  echo "BWA indexes a fasta reference sequence"
  echo "The resulting indicies will be placed in the same directory as the reference sequence"
  echo
  echo "  -r: Path to fasta reference sequence"
  echo
  exit 1
}

while getopts "r:" opt; do
  case "${opt}" in
    r )
        ref_file=${OPTARG}
      ;;
    \? )
        display_help
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$ref_file" ]]; then
    display_help
    exit 1
fi

ref_dir=$(dirname ${ref_file})
ref_name=$(basename ${ref_file})

cp -vL ${ref_dir}/${ref_name} $TMPDIR
cd $TMPDIR
bwa index ${ref_name}
samtools faidx ${ref_name} > ${ref_name}.fai
samtools dict -o ${ref_name}.dict ${ref_name}
rm $TMPDIR/${ref_name}
cd -
cp -v $TMPDIR/* ${ref_dir}/

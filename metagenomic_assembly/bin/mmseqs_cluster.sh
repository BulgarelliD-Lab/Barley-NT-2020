#!/bin/env bash

#$ -cwd
#$ -pe smp 48
#$ -jc long
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID

# NB This script uses $TMPDIR1 which is mapped to a ramdisk for performance
# If unavailble change to $TMPDIR which uses on-disk temporary space

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -i [input_fasta] -o [output_dir]"
  echo
  echo "Carried out clustering of protein predictions using mmseqs2 'easy-cluster' method"
  echo
  echo "  -i: input fasta file."
  echo "  -o: output directory"
  echo
  exit 1
}
while getopts "i:o:" opt; do
  case "${opt}" in
    i )
        in_file=${OPTARG}
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

if [[ -z "$in_file" ]] || [[ -z "$out_dir" ]]; then
    display_help
    exit 1
fi

echo "Hostname=$HOSTNAME"
echo "Input file=${in_file}"
echo "Output directory=${out_dir}"

if [ ! -d "${out_dir}" ]; then
	mkdir -pv ${out_dir}
fi

mmseqs easy-cluster ${in_file} ${out_dir}/mmseqs $TMPDIR1 --threads 48


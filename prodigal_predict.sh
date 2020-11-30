#!/bin/env bash

#$ -cwd 
#$ -V
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail
display_help() {
  echo "Usage: $(basename $0) -i [input_dir] -o [output_dir]"
  echo
  echo "Carries out prediction of protein sequences from an assembly using prodigal"
  echo
  echo "  -i: input assembly directory. Should contain a separate subdirectory for each"
  echo "	  metaspades assembly, containing a 'contigs.fasta' file."
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
    display_help
    exit 1
fi

if [[ -z "$SGE_TASK_ID" ]] || [[ "$SGE_TASK_ID" == 'undefined' ]]; then
    echo "Error: Script must be submitted as an array job (qsub -t 1-n)"
    exit 1
fi

mkdir -p ${out_dir}

readarray -t samples < <(ls -1 ${in_dir})
sample=${samples[$SGE_TASK_ID-1]}

echo "Hostname=$HOSTNAME"
echo "Sample=${sample}"
echo "Input directory=${in_dir}"
echo "Output directory=${out_dir}"

cp -v ${in_dir}/${sample}/contigs.fasta $TMPDIR

prodigal -i $TMPDIR/contigs.fasta -o $TMPDIR/${sample}.proteins.gff3 -d $TMPDIR/${sample}.nucl.fa -a $TMPDIR/${sample}.proteins.fa -p meta -f gff

cp -v $TMPDIR/${sample}.proteins.* ${out_dir}
cp -v $TMPDIR/${sample}.nucl.* ${out_dir}


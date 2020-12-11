#!/bin/env bash

#$ -cwd
#$ -pe smp 24
#$ -mods l_hard mfree 96G
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

set -eo pipefail

usage="Usage: $0 -i /path/to/input/fastq/directory -o /path/to/output/directory -d path/to/kraken_db" 
display_help() {
  echo "Usage: $(basename $0) -i fastq_directory -o output_directory -d kraken_database"
  echo 
  echo  "Runs kraken2 against paired fastq reads, where read filenames start with the sample followed by a '_'."
  echo
  echo  "This script is designed to be run as an array job on an SGE cluster"
  echo  "which has local disk storage available on execution nodes accessed via $TMPDIR"
  echo
  echo "  -f: input fastq directory. Files should be paired with pair indicated by _[12]"
  echo "      and gzip compressed"
  echo "  -o: output directory"
  echo "  -d: kraken database"
  echo
  exit 1
}

while getopts "i:o:d:" opt; do
  case "${opt}" in
    i )
    	in_dir=${OPTARG}
    	;;
    o )
      out_dir=${OPTARG}
      ;;
    d )
      db=${OPTARG}
      ;;
    \? )
      display_help
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$in_dir" ]] || [[ -z "$out_dir" ]] || [[ -z "$db" ]]; then
    display_help
    exit 1
fi

if [[ ! -e "${in_dir}" ]]; then
	echo "Specificed input directory (${indir}) does not exist"
	exit
fi

if [[ ! -e "${db}" ]]; then
	echo "Specificed kraken database (${db}) does not exist"
	exit
fi

if [[ -z "${SGE_TASK_ID}" ]]; then
  echo "This script should be executed as an SGE array job"
  exit
fi

mkdir -p ${out_dir}

readarray -t samples < <(ls $in_dir/*gz|cut -f1 -d_|uniq)
sample=$(basename ${samples[$SGE_TASK_ID-1]})
readarray -t files < <(ls -1 $in_dir/${sample}*|xargs -i basename {})

mkdir $TMPDIR/kraken_db
cp -v ${db}/*k2d $TMPDIR/kraken_db
cp -v ${in_dir}/${sample}* $TMPDIR

kraken2 --db $TMPDIR/kraken_db --threads 24 --use-names --report ${out_dir}/${sample}.report.txt --output ${out_dir}/${sample}.txt --gzip-compressed --paired ${TMPDIR}/${files[0]} ${TMPDIR}/${files[1]}



#!/bin/env bash

#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -cwd 
#$ -pe smp 48
#$ -jc long
#$ -mods l_hard h_vmem 12G
#$ -adds l_hard ramdisk 110G

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -i [input_fasta] -o [output_dir] -e eggnog_data_dir"
  echo
  echo "Functional annotation of preticted protein clusters using eggnog_mapper."
  echo
  echo "Eggnog data (eggnog.db.gz and eggnog_proteins.dmnd) must be provided in a"
  echo "directory which is passed using the '-e' argument"
  echo
  echo "This script uses $TMPDIR1 which is assigned to a ram-disk for performance"
  echo "purposes, and $TMPDIR for disk-based temporary data"
  echo
  echo "  -i: input fasta file."
  echo "  -o: output directory"
  echo "  -e: directory containing eggnog data"
  echo
  exit 1
}

while getopts "i:e:o:" opt; do
  case "${opt}" in
    i ) 
		input=${OPTARG}
      ;;
	o )
		out_dir=${OPTARG}
		;;
	e )
		egg_data=${OPTARG}
	;;
    \? ) 
		echo ${usage}
      ;;
	: )
		echo "Invalid option: $OPTARG requires an argument" 1>&2
	;;
  esac
done

if [[ -z "${input}" ]] || [[ -z "${out_dir}" ]] || [[ -z "${egg_data}" ]]; then
	display_help
	exit 1
fi

if [[ ! -d ${out_dir} ]]; then
	mkdir -pv ${out_dir}
fi

echo hostname=${HOSTNAME}
echo input=${input}
echo out_dir=${out_dir}
echo egg_data=${egg_data}

if [[ ! -f "${egg_data}/eggnog.db.gz" ]] || [[ ! -f "${egg_data}/eggnog_proteins.dmnd" ]]; then
	echo "Required eggnog databases not found."
	echo
	echo "Please ensure eggnog.db.gz and eggnog_proteins.dmnd"
	echo "are present in ${egg_data}"
	exit 1
fi

gunzip -c ${egg_data}/eggnog.db.gz > ${TMPDIR1}/eggnog.db
cp ${egg_data}/eggnog_proteins.dmnd ${TMPDIR1}/

emapper.py -i ${input} -o eggnog --output_dir ${out_dir} --scratch_dir ${TMPDIR} --data_dir ${TMPDIR1} \
	--database bact -m diamond --cpu 48 --go_evidence non-electronic --target_orthologs all \
	--seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --override


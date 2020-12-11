#!/bin/env bash

#$ -cwd
#$ -V
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -pe smp 8
#$ -mods l_hard m_mem_free 16G
#$ -adds l_hard h_vmem 16G

set -eo pipefail

display_help() {
  echo "Usage: $(basename $0) -i [input_fasta] -o [output_dir]"
  echo
  echo "Separates alignments into separate BAM files containing"
  echo "properly paired alignments, singletons and alignments"
  echo "split between separate reference sequences"
  echo
  echo "  -d: directory containing BAM alignments"
  echo
  exit 1
}

while getopts "d:" opt; do
  case "${opt}" in
    d )
        dir=${OPTARG}
      ;;
    \? )
        echo ${usage}
      ;;
    : )
        echo "Invalid option: $OPTARG requires an argument" 1>&2
    ;;
  esac
done

if [[ -z "$dir" ]]; then
	display_help
	exit 1
fi

if [[ -z "$SGE_TASK_ID" ]] || [[ "$SGE_TASK_ID" == 'undefined' ]]; then
	echo "This script must be submitted as an array job"
	exit 1
fi

readarray -t files < <(ls -1 ${dir}/*bam)
bamfile=$(basename ${files[$SGE_TASK_ID-1]})
sample=$(basename ${files[$SGE_TASK_ID-1]%%'.'*})

echo hostname=${HOSTNAME}
echo dir=${dir}
echo sample=${sample}
echo bamfile=${bamfile}

cp -v ${dir}/${bamfile} $TMPDIR

samtools index -@ 8 $TMPDIR/${bamfile}
# Proper pairs
samtools view -@ 8 -q 30 -f 2 -F 256 -b -o $TMPDIR/${sample}.proper.bam $TMPDIR/${bamfile}
# Diff cDNA
samtools view -H $TMPDIR/${bamfile} > $TMPDIR/${sample}.diff_cdna.sam
samtools view -@ 8 -q 30 $TMPDIR/${bamfile} | awk '($3!=$7 && $7!="=")' >> $TMPDIR/${sample}.diff_cdna.sam
samtools view -@ 8 -b -o $TMPDIR/${sample}.diff_cdna.bam $TMPDIR/${sample}.diff_cdna.sam
# Singletons
samtools view -@ 8 -q 30 -f 8 -F 256 -o $TMPDIR/${sample}.singletons.bam $TMPDIR/${bamfile}

rm $TMPDIR/${bamfile}
rm $TMPDIR/${sample}.diff_cdna.sam
for type in proper diff_cdna singletons; do
  samtools index -@ 8 $TMPDIR/${sample}.${type}.bam
  samtools flagstat -@ 8 $TMPDIR/${sample}.${type}.bam > $TMPDIR/${sample}.${type}.flagstat
  samtools idxstats -@ 8 $TMPDIR/${sample}.${type}.bam > $TMPDIR/${sample}.${type}.idxstats
done
cp -v $TMPDIR/* ${dir}



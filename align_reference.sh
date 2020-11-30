#!/bin/env bash

#$ -cwd 
#$ -V
#$ -pe smp 24
#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID.$TASK_ID

# Carries out alignment against a reference genome, sorting of outputs
# and marking duplicate reads

# Resulting outputs are written to an 'alignments' directory of the cwd

# Requires BWA and picard tools on the path

set -eo pipefail
sample=""

usage="Usage: $0 -d [reference_dir] -r [reference_name] -f [fastq_dir] -o [output_dir] [-s sample]"
display_help() {
  echo "Usage: $(basename $0) -d reference_dir -r reference_name -f fastq_dir -o output_dir [-s sample]"
  echo
  echo "Carries out alignment against a reference genome, sorting of outputs and marking duplicate reads"
  echo "Fastq reads should be in a separate subdirectory of the 'fastq_dir', named with the sample ID"
  echo 
  echo "If submitted as an array job, each task operates on a separate sample"
  echo "Optionally, a single sample can be specified with the '-s' argument and the job submitted as 'normal' job"
  echo
  echo "  -f: input fastq directory. Files should be paired and suffixed with _[12].fq.gz"
  echo "  -o: output directory"
  echo "  -d: directory containing BWA indices"
  echo "  -r: name of reference index"
  echo "  -s: sample id [optional]"
  echo
  exit 1
}

while getopts "d:r:f:o:s:" opt; do
  case "${opt}" in
    d )
        ref_dir=${OPTARG}
      ;;
    r )
        ref_name=${OPTARG}
      ;;
    s )
        sample=${OPTARG}
      ;;
	  f ) 
		    fastq_dir=${OPTARG}
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

if [[ -z "$ref_dir" ]] || [[ -z "$ref_name" ]] || [[ -z "$fastq_dir" ]] || [[ -z "$out_dir" ]]; then
    display_help
    exit 1
fi

if [[ "$SGE_TASK_ID" == 'undefined' ]] && [[ -z ${sample} ]]; then
    echo "A sample must be specified with '-s' if submitting as a non-array job"
    exit 1
fi

mkdir -p ${out_dir} 

readarray -t samples < <(ls -1 ${fastq_dir}|cut -f1 -d_|sort -u)
if [[ -z "$sample" ]]; then
  sample=${samples[$SGE_TASK_ID-1]}
fi
echo
echo "Hostname=${HOSTNAME}"
echo "Sample=${sample}"
echo "ref_dir=${ref_dir}"
echo "ref_name=${ref_name}"
echo "fastq_dir=${fastq_dir}"
echo "out_dir=${out_dir}"
echo

cp -Lrv ${ref_dir}/${ref_name}* $TMPDIR
cp -Lv ${fastq_dir}/${sample}* $TMPDIR
r1=$(ls $TMPDIR/${sample}*1.fq.gz)
r2=$(ls $TMPDIR/${sample}*2.fq.gz)

echo "r1=$r1"
echo "r2=$r2"
echo

bwa mem -t 24 -M $TMPDIR/${ref_name} ${r1} ${r2} \
  |samtools view -b -o $TMPDIR/${sample}.bam

samtools sort -@ 24 -o $TMPDIR/${sample}.sorted.bam $TMPDIR/${sample}.bam
samtools index -c -m 14 -@ 24 $TMPDIR/${sample}.sorted.bam

picard MarkDuplicates I=$TMPDIR/${sample}.sorted.bam O=$TMPDIR/${sample}.sorted.dup.bam \
  M=$TMPDIR/${sample}.duplicates.txt VALIDATION_STRINGENCY=LENIENT

samtools index -c -m 14 -@ 24 $TMPDIR/${sample}.sorted.dup.bam
samtools flagstat -@ 8 $TMPDIR/${sample}.sorted.dup.bam > $TMPDIR/${sample}.sorted.dup.flagstat
samtools idxstats -@ 8 $TMPDIR/${sample}.sorted.dup.bam > $TMPDIR/${sample}.sorted.dup.idxstats

cp -v $TMPDIR/*sorted.dup* ${out_dir}

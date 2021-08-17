#!/user/bin/bash

bam_files=$1

conda activate winlow

set -- $(ls -1 "$bam_files"/*.bam | grep -v 'tagged' | sort -t - -k2n -k1V)
[[ -e $1 || -L $1 ]] || { echo "No .bam files found in $bam_files" >&2; exit 1; }

while (( $# > 1 ))
do

  stem=$(basename "${1}" )
  output_base=$(echo $stem | cut -d '.' -f 1)

  echo "Running winlow for $output_base"

  echo "python winlow.py -r 3R:0-3000000 -w 100000 -s 50 --tumour $1 --normal $2 --debug"

  python winlow.py -r 3R:0-3000000 -w 100000 -s 50 --tumour $1 --normal $2 --debug

  # echo "bash $script $output_base $n_base $1 $2 $group $out_dir $sex"
  # bash $script $output_base $n_base $1 $2 $group $out_dir $sex >> $submit_log
  # config=$out_dir/${output_base}.config
  #
  # qsub -v VAR1=$config,VAR2=$out_dir,VAR3=$output_base -o $log/${output_base}.runlog -j oe -N ${output_base}.FREEC $pbs_dir/run_freec.pbs
  # echo "DONE" >> $submit_log

  shift 2 || break

done

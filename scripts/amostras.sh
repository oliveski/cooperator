#!/bin/bash

r=0.0
i=1001	# primeira semente
samples=100

cp job.sh START.sh
echo "echo \"It is now $(date)\" > START" >> START.sh
sbatch START.sh
rm START.sh

mkdir JOBS
while (( $(bc <<<"$r <= 3") )); do
	count=0
	while [ $count -lt $samples ]; do
		run_file=run$i.sh
		cp job.sh $run_file
		echo "echo \"eu sou o $r $i\"" >> $run_file
		echo "time ./exe 50 $r $i 0.8 0.5" >> $run_file
		sbatch $run_file
		mv $run_file JOBS
		i=$(($i + 2))
		count=$(($count + 1))
	done
	r=$(bc <<< "$r + 0.1")
done

cp job.sh END.sh
echo "echo \"It is now $(date)\" > END" >> END.sh
sbatch END.sh
rm END.sh

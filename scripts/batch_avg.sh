#!/bin/bash

# Takes averages over a batch of samples

### script_dir is where the awk script is
### data_dir is where the data files are

script_dir=

data_dir=

samples=
count=0
for file in $data_dir/*.dat; do
	echo "$file $($script_dir/media.awk $file)" >> $data_dir/tmp.file
	echo $file
	count=$(($count + 1))
	if [[ count -eq samples ]]; then
		echo "$file $($script_dir/media.awk $data_dir/tmp.file)" >> $data_dir/medias.dsf
		rm $data_dir/tmp.file
		count=0
		echo ""
	fi
done


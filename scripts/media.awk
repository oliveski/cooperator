#!/bin/awk -f
# Change the col to be averaged over
# with col and change the number of
# header lines on num_head
BEGIN{
	sum = 0; summ = 0;
	num_head = 0;
	col = 2;
}
{
	if(NR != num_head){
		arr[NR] = $(col); sum += $(col);
	}
}
END{
	if(NR == 1){
		print $(col) " " 0;
		exit 0;
	}
	media = sum / (NR - num_head);
	for(i = num_head + 1; i <= NR; i++){
		summ = summ + (arr[i] - media)^2;
	};
	summ = sqrt(summ/(NR - 1 - num_head));
	print media " " summ;
}

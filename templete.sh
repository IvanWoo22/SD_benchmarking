for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 0 "\t" $1}' >time1-0.txt
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time1-1.txt
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (i = (a+1); i <= NF; i++) {
					if ($i == 0) {
							count0++;
					} else if ($i == 1) {
							count1++;
					} else if ($i == 2) {
							count2++;
					} else if ($i == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time1-2.txt
		for k in {2..5}; do
			cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
				| awk -va="$k" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (k = (a+1); k <= NF; k++) {
						if ($k == 0) {
								count0++;
						} else if ($k == 1) {
								count1++;
						} else if ($k == 2) {
								count2++;
						} else if ($k == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$k" '$a==3{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time"${k}"-1.txt
			cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
				| awk -va="$k" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (k = (a+1); k <= NF; k++) {
						if ($k == 0) {
								count0++;
						} else if ($k == 1) {
								count1++;
						} else if ($k == 2) {
								count2++;
						} else if ($k == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$k" '$a==3{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time"${k}"-2.txt
			cat time"${k}"-{1,2}.txt >${i}_${j}/time"${k}".txt
		done
		cat time1-{1,2}.txt >${i}_${j}/time1.txt
		cut -f 3,4,5,6,7 ${i}_${j}/time-point.tsv | sort | uniq -c \
			| awk '$NF==3{print $2 $3 $4 $5 $6 $7 "\t" 0 "\t" $1}' >${i}_${j}/time6.txt
		rm time*.txt
		for k in {1..6}; do
			perl ../fetch_timespot.pl ${i}_${j}/time"${k}".txt ${i}_${j}/time-point.tsv \
				>${i}_${j}/time-point."${k}".tsv
			awk '{print $1 "\t" $2 "\t" $3 "\t" NR ",1\n" $4 "\t" $5 "\t" $6 "\t" NR ",2"}' \
				${i}_${j}/time-point."${k}".tsv \
				| sort -k1,1 -k2,2n >${i}_${j}/time-point."${k}".sort.bed
			perl ../arg_meth_link.pl \
				../../MASED/Memory/AT.beta.1.tsv \
				${i}_${j}/time-point."${k}".sort.bed \
				${i}_${j}/time-point."${k}".beta.bed
		done
	done
done

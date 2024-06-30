perl ../fetch_time_point.pl \
	{1}_{2}/time{3}.txt {1}_{2}/time-point.tsv \
	>{1}_{2}/time-point.{3}.tsv
awk '$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2}
$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0}' \
	{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n \
	>{1}_{2}/time-point.{3}.sort.bed
awk '$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2" $7 "\t" $9}
$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0" $7 "\t" $9}' \
	{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.evolution.bed
parallel -j 24 '
	perl ../fetch_time_point.pl {1}_{2}/time{3}.txt {1}_{2}/time-point.tsv >{1}_{2}/time-point.{3}.tsv
	awk '\''$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2"}
					$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0"}'\'' \
		{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.sort.bed
	awk '\''$8==1{print $1 "\t" $2 "\t" $3 "\t" NR "\t1\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t2\t" $7 "\t" $9}
					$8==0{print $1 "\t" $2 "\t" $3 "\t" NR "\t0\t" $7 "\t" $9 "\n" $4 "\t" $5 "\t" $6 "\t" NR "\t0\t" $7 "\t" $9}'\'' \
		{1}_{2}/time-point.{3}.tsv | sort -k1,1 -k2,2n >{1}_{2}/time-point.{3}.evolution.bed
' ::: lastz biser ::: TAIR10_rmmasked{1..5} ::: {1..6}

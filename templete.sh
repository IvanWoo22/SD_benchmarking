for k in {1..6}; do
	closestBed -d -a \
		../sd_meth/biser_TAIR10_rmmasked5/time-point."${k}".evolution.bed \
		-b ../sd_meth/TE.bed -t all \
		| awk '{
			if ($NF == 0) {
        if ($9 <= $2) {
            start = 0
        } else {
            start = $9 - $2
        }
        if ($10 >= $3) {
            end = $3 - $2
        } else {
            end = $10 - $2
        }
        printf "Chr%s(+):%d-%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, start, end, $4, $5, $6, $7, $11
			}
	}' >TE.list
	perl ../region_blank_pick.pl \
		-f biser_TAIR10_rmmasked5/time-point."${k}".muscle.fa -r TE.list \
		>biser_TAIR10_rmmasked5/time-point."${k}".muscle.TE.tsv
done

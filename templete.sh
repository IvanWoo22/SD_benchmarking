for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		for k in {1..6}; do
			awk -va=${i} -vb=${j} -vc="${k}" '{sum+=$3};END{print b " " a "\t" c "\t" sum}' ${i}_${j}/time"${k}".txt >>links2_timedistribution.tsv
		done
		awk -F"\t" -va=${i} -vb=${j} '$1==b" "a{sum+=$3};END{print b " " a "\t" sum}' links2_timedistribution.tsv >>links2_timedistributionallcount.tsv
		wc -l <${i}_${j}/links.2copy.tsv >>links2_allcount.tsv
	done
done

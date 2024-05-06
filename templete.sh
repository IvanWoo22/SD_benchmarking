rm time-point.prom.beta.all.bed
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.prom.beta.bed
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.prom.beta.bed >>${i}_${j}/time-point.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.prom.beta.bed >>time-point.prom.beta.all.bed
		done
	done
done

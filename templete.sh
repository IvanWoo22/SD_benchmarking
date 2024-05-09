parallel -j 20 "
	perl ../promotor_intsec.pl \
		../struct/pseudogene.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.pseudo.prom.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.pseudo.prom.bed \
		{1}_{2}/time-point.{3}.pseudo.prom.beta.bed 6
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.pseudo.prom.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.pseudo.prom.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pseudo.prom.beta.bed \
				>>${i}_${j}/time-point.pseudo.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pseudo.prom.beta.bed \
				>>time-point.pseudo.prom.beta.all.bed
		done
	done
done

parallel -j 20 "
	perl ../promotor_intsec.pl \
		../struct/protein_coding.promoter.bed {1}_{2}/time-point.{3}.sort.bed >{1}_{2}/time-point.{3}.pc.prom.bed
	perl ../arg_meth_link_neo.pl \
		../../MASED/Memory/AT.beta.1.tsv \
		{1}_{2}/time-point.{3}.pc.prom.bed \
		{1}_{2}/time-point.{3}.pc.prom.beta.bed 6
" ::: lastz biser ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ::: {1..6}

rm time-point.pc.prom.beta.all.bed 2>/dev/null
for i in lastz biser; do
	for j in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
		rm ${i}_${j}/time-point.pc.prom.beta.bed 2>/dev/null
		for k in {1..6}; do
			awk -va="${k}" '{print $1 "\t" $2 "\t" $3 "\t" a "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pc.prom.beta.bed \
				>>${i}_${j}/time-point.pc.prom.beta.bed
			awk -va=${i} -vb=${j} -vc="${k}" '{print b "\t" a "\t" $1 "\t" $2 "\t" $3 "\t" c "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $6}' \
				${i}_${j}/time-point.${k}.pc.prom.beta.bed \
				>>time-point.pc.prom.beta.all.bed
		done
	done
done

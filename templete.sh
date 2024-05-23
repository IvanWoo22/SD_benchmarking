for i in TAIR10 ColCEN; do
	mkdir ${i}_rm && cd ${i}_rm || exit
	faops split-name ../${i}_unmasked.fa .
	parallel -j 5 '
		RepeatMasker -species arabidopsis -pa 5 -s -xsmall -e ncbi -dir . Chr{}.fa
		RM2Bed.py -d . Chr{}.fa.out
		trf Chr{}.fa 2 5 7 80 10 40 500 -l 10 >Chr{}.fa.dat
		dustmasker -in Chr{}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr{}.fa_dust.bed
	' ::: {1..5}
	cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
	cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
	python ../../trf_merge.py Chr{1..5}.fa.dat trf.out.bed

	cat trf.out.bed dust.out.bed repeatmasker.out.bed \
		| cut -f 1-3 >>tmp.msk.bed
	cut -f 1-3 tmp.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -i - \
		| awk '$3-$2 > 2000 {print $0}' \
		| bedtools merge -d 100 -i - >tmp1.msk.bed
	cut -f 1-3 tmp.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -i - \
		| awk '$3-$2 <= 2000 {print $0}' \
		| bedtools merge -d 10 -i - >tmp2.msk.bed
	cut -f 1-3 tmp.msk.bed tmp1.msk.bed tmp2.msk.bed \
		| bedtools sort -i - \
		| bedtools merge -i - \
		| seqtk seq -l 50 -M /dev/stdin ../${i}_unmasked.fa >../${i}_rmmasked.fa
	cd ..
done

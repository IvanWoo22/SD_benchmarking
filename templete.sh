mkdir "biser" && cd "biser" || exit
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	biser -t 20 -o ${PREFIX}.out ../reference/${PREFIX}.fa
	awk '{print $1"("$9"):"$2"-"$3"\t"$4"("$10"):"$5"-"$6}' ${PREFIX}.out \
		| linkr sort stdin \
		| linkr clean stdin -o ${PREFIX}.links.sort.clean.tsv
	linkr connect ${PREFIX}.links.sort.clean.tsv -r 0.05 \
		| linkr filter stdin -r 0.05 -o ${PREFIX}.links.filter.tsv
	linkr filter ${PREFIX}.links.filter.tsv -n 2 -o stdout >${PREFIX}.links2.tsv
done
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	rgr merge ${PREFIX}.links.sort.clean.tsv -c 0.95 -o ${PREFIX}.links.merge.tsv
	linkr clean ${PREFIX}.links.sort.clean.tsv -r ${PREFIX}.links.merge.tsv --bundle 500 -o ${PREFIX}.links.clean.tsv
done
cd ..

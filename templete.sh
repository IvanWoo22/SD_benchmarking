for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	linkr connect ${PREFIX}.links.sort.clean.tsv -r 0.05 \
		| linkr filter stdin -r 0.05 -o ${PREFIX}.links.filter.tsv
done
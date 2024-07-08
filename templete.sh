perl ../bed_intsec.pl \
	--diff --ref_bed ./TE.bed \
	--in_bed {1}_{2}/time-point.{3}.sort.bed \
	--output {1}_{2}/time-point.{3}.nonTE.bed
perl ../arg_meth_link_neo.pl \
	../../MASED/Memory/AT.beta.1.tsv \
	{1}_{2}/time-point.{3}.nonTE.bed \
	{1}_{2}/time-point.{3}.nonTE.beta.bed 6

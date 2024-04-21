for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	egaz prepseq ../reference/${PREFIX}.fa -o ${PREFIX}/
	cd ${PREFIX} || exit
	perl -nla -F"\t" -e "print for @F" <links.filter.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.filter.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv \
			| perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " \
				>>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	mv cover.yml Atha.cover.yml
	mv copy.yml Atha.copy.yml
	mv cover.yml.csv Atha.cover.csv
	mv copy.csv Atha.copy.csv
	mv links.filter.tsv Atha.links.tsv
	rm -rf Chr*.fa chr* ./*.temp.yml ./*.links.merge.tsv ./*.links.sort.clean.tsv
	cd ..
done

for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_rmmasked; do
	cd ${PREFIX}/Processing/Atha || exit
	egaz blastn axt.gl.fasta genome.fa -o axt.bg.blast --parallel 8
	egaz blastmatch axt.bg.blast -c 0.95 -o axt.bg.region --parallel 8
	samtools faidx genome.fa -r axt.bg.region --continue \
		| perl -p -e "/^>/ and s/:/(+):/" >axt.bg.fasta
	cat axt.gl.fasta axt.bg.fasta | faops filter -u stdin stdout \
		| faops filter -n 250 stdin stdout >axt.all.fasta
	egaz blastn axt.all.fasta axt.all.fasta -o axt.all.blast --parallel 8
	egaz blastlink axt.all.blast -c 0.95 -o links.blast.tsv --parallel 8

	linkr sort links.lastz.tsv links.blast.tsv -o links.sort.tsv
	linkr clean links.sort.tsv -o links.sort.clean.tsv

	rgr merge links.sort.clean.tsv -c 0.95 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 500 -o links.clean.tsv

	linkr connect links.clean.tsv -r 0.9 -o links.connect.tsv
	linkr filter links.connect.tsv -r 0.8 -o links.filter.tsv

	fasops create links.filter.tsv -g genome.fa -o multi.temp.fas
	fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o multi.refine.fas
	fasops links multi.refine.fas -o stdout | linkr sort stdin -o stdout | linkr filter stdin -n 2-50 -o links.refine.tsv
	fasops links multi.refine.fas -o stdout --best | linkr sort stdin -o links.best.tsv
	fasops create links.best.tsv -g genome.fa --name Atha -o pair.temp.fas
	fasops refine pair.temp.fas --msa mafft -p 16 -o pair.refine.fas

	perl -nla -F"\t" -e "print for @F" <links.refine.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.refine.tsv -n ${n} -o stdout \
			>links.copy${n}.tsv
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
	cp cover.yml ../../Atha.cover.yml
	cp copy.yml ../../Atha.copy.yml
	mv cover.yml.csv ../../Atha.cover.csv
	mv copy.csv ../../Atha.copy.csv
	cp links.refine.tsv ../../Atha.links.tsv
	mv multi.refine.fas ../../Atha.multi.fas
	mv pair.refine.fas ../../Atha.pair.fas
	cd ../..
	rm -rf Processing Pairwise Results Chr*.fa chr* *.sh
	cd ..
done

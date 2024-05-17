mkdir "preref" && cd "preref" || exit
cp ../data/TAIR10_chr_all.fas.gz .
cp ../data/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz .
cp ../data/Col-CEN_v1.2.fasta.gz .
pigz -dc TAIR10_chr_all.fas.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >TAIR10_unmasked.fa
rm Chr*.fa
pigz -dc Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz | sed 's/>/>Chr/g' | faops split-name stdin .
rm ChrPt.fa ChrMt.fa
cat Chr*.fa >TAIR10_E58masked.fa
rm Chr*.fa
pigz -dc Col-CEN_v1.2.fasta.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >ColCEN_unmasked.fa
rm Chr*.fa
rm TAIR10_chr_all.fas.gz Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz Col-CEN_v1.2.fasta.gz

parallel -j 2 '
	mkdir {}_rm && cd {}_rm || exit
	faops split-name ../{}_unmasked.fa .
	for j in 1 2 3 4 5; do
		RepeatMasker -species arabidopsis -pa 12 -s -xsmall -e ncbi -dir . Chr${j}.fa
		RM2Bed.py -d . Chr${j}.fa.out
		trf Chr${j}.fa 2 7 7 80 10 50 15 -l 25 -h -ngs >Chr${j}.fa.dat
		dustmasker -in Chr${j}.fa -outfmt acclist -out - | sed '\''s/^>//'\'' >Chr${j}.fa_dust.bed
	done
	cat Chr*.fa_rm.bed | bedtools sort -i - >repeatmasker.out.bed
	cat Chr*.fa_dust.bed | bedtools sort -i - >dust.out.bed
	python ../../trf_merge.py Chr{1..5}.fa.dat trf.out.bed

	cat trf.out.bed dust.out.bed repeatmasker.out.bed | cut -f 1-3 >>tmp.msk.bed
	cut -f 1-3 tmp.msk.bed | bedtools sort -i - | bedtools merge -i - | awk '\''$3-$2 > 2000 {print $0}'\'' | bedtools merge -d 100 -i - >tmp2.msk.bed
	cut -f 1-3 tmp.msk.bed tmp2.msk.bed | bedtools sort -i - | bedtools merge -i - | seqtk seq -l 50 -M /dev/stdin ../{}_unmasked.fa >../{}_rmmasked.fa
' ::: TAIR10 ColCEN

echo ",chr length,masked size,masked coverage" >masked_cover.csv
mkdir "smcover" && cd "smcover" || exit
egaz prepseq ../TAIR10_unmasked.fa -o .
for i in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked; do
	{
		echo -ne "${i},"
		faops masked ../${i}.fa | spanr cover stdin \
			| spanr stat --all chr.sizes stdin \
			| tail -1
	} >>../masked_cover.csv
done
rm chr.* Chr*
egaz prepseq ../ColCEN_unmasked.fa -o .
{
	awk '{sum+=$2};END{print "ColCEN_unmasked," sum ",0,0.0000"}' chr.sizes
	echo -ne "ColCEN_rmmasked,"
	faops masked ../ColCEN_rmmasked.fa | spanr cover stdin \
		| spanr stat --all chr.sizes stdin \
		| tail -1
} >>../masked_cover.csv
cd .. && rm -rf "smcover"
cd ..

mkdir "biser" && cd "biser" || exit
parallel -j 6 '
	mkdir {} && cd {} || exit
	biser -t 4 -o biser_out ../../preref/{}.fa
	awk '\''{print $1"("$9"):"$2"-"$3"\t"$4"("$10"):"$5"-"$6}'\'' biser_out | linkr sort stdin | linkr clean stdin -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.95 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 500 -o links.clean.tsv
	linkr connect links.clean.tsv -r 0.05 | linkr filter stdin -r 0.05 -o links.tsv
	singularity run $HOME/egaz_master.sif egaz prepseq ../../preref/{}.fa -o .
	perl -nla -F"\t" -e "print for @F" <links.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv | perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " >>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	rm -rf Chr*.fa chr* ./*.temp.yml ./links.merge.tsv ./links.sort.clean.tsv
' ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked
cd ..

mkdir "lastz" && cd "lastz" || exit
echo 'strain,strain_id,species,species_id,genus,genus_id,family,family_id,order,order_id
Atha,3702,"Arabidopsis thaliana",3702,Arabidopsis,3701,Brassicaceae,3700,Brassicales,3699' >ensembl_taxon.csv
parallel -j 6 '
	mkdir {}
	singularity run $HOME/egaz_master.sif egaz prepseq ../preref/{}.fa -o {}/
	singularity run $HOME/egaz_master.sif egaz template {}/ --self -o {}/ --taxon ./ensembl_taxon.csv --circos --parallel 16 -v
	cd {} || exit
	mkdir -p Pairwise
	singularity run $HOME/egaz_master.sif egaz lastz --isself --set set01 -C 0 --parallel 4 --verbose . . -o Pairwise/vsSelf
	singularity run $HOME/egaz_master.sif egaz lpcnam --parallel 4 --verbose . . Pairwise/vsSelf
	mkdir -p Results/vsSelf
	mkdir -p Processing/vsSelf
	ln -s "$(pwd)"/chr.fasta Processing/vsSelf/genome.fa
	cp -f "$(pwd)"/chr.sizes Processing/vsSelf/chr.sizes
	cd Processing/vsSelf || exit
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops axt2fas ../../Pairwise/vsSelf/axtNet/*.axt.gz -l 1000 -s chr.sizes -o stdout >axt.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops separate axt.fas -o . --nodash -s .sep.fasta
	singularity run $HOME/egaz_master.sif egaz exactmatch target.sep.fasta genome.fa --length 500 --discard 50 -o replace.target.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops replace axt.fas replace.target.tsv -o axt.target.fas
	singularity run $HOME/egaz_master.sif egaz exactmatch query.sep.fasta genome.fa --length 500 --discard 50 -o replace.query.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops covers axt.correct.fas -o axt.correct.yml
	spanr split axt.correct.yml -s .temp.yml -o .
	spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
	spanr stat chr.sizes axt.union.yml -o union.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops links axt.correct.fas -o stdout | perl -nl -e "s/(target|query)\.//g; print;" >links.lastz.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops separate axt.correct.fas --nodash --rc -o stdout | perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" | faops filter -u stdin stdout | faops filter -n 250 stdin stdout >axt.gl.fasta
	singularity run $HOME/egaz_master.sif egaz blastn axt.gl.fasta genome.fa -o axt.bg.blast --parallel 8
	singularity run $HOME/egaz_master.sif egaz blastmatch axt.bg.blast -c 0.95 -o axt.bg.region --parallel 8
	samtools faidx genome.fa -r axt.bg.region --continue | perl -p -e "/^>/ and s/:/(+):/" >axt.bg.fasta
	cat axt.gl.fasta axt.bg.fasta | faops filter -u stdin stdout | faops filter -n 250 stdin stdout >axt.all.fasta
	singularity run $HOME/egaz_master.sif egaz blastn axt.all.fasta axt.all.fasta -o axt.all.blast --parallel 8
	singularity run $HOME/egaz_master.sif egaz blastlink axt.all.blast -c 0.95 -o links.blast.tsv --parallel 8
	linkr sort links.lastz.tsv links.blast.tsv -o links.sort.tsv
	linkr clean links.sort.tsv -o links.sort.clean.tsv
	rgr merge links.sort.clean.tsv -c 0.95 -o links.merge.tsv
	linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 500 -o links.clean.tsv
	linkr connect links.clean.tsv -r 0.9 -o links.connect.tsv
	linkr filter links.connect.tsv -r 0.8 -o links.filter.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops create links.filter.tsv -g genome.fa -o multi.temp.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o multi.refine.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops links multi.refine.fas -o stdout | linkr sort stdin -o stdout | linkr filter stdin -n 2-50 -o links.refine.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops links multi.refine.fas -o stdout --best | linkr sort stdin -o links.best.tsv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops create links.best.tsv -g genome.fa --name Atha -o pair.temp.fas
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops refine pair.temp.fas --msa mafft -p 16 -o pair.refine.fas
	perl -nla -F"\t" -e "print for @F" <links.refine.tsv | spanr cover stdin -o cover.yml
	echo "key,count" >links.count.csv
	for n in 2 3 4-50; do
		linkr filter links.refine.tsv -n ${n} -o links.copy${n}.tsv
		perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
		wc -l links.copy${n}.tsv | perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " >>links.count.csv
		rm links.copy${n}.tsv
	done
	spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
	spanr stat chr.sizes copy.yml --all -o links.copy.csv
	/home/linuxbrew/.linuxbrew/Cellar/perl/5.38.2_1/bin/fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
	spanr stat chr.sizes cover.yml -o cover.yml.csv
	cp cover.yml ../../cover.yml
	cp copy.yml ../../copy.yml
	mv cover.yml.csv ../../cover.csv
	mv copy.csv ../../copy.csv
	cp links.refine.tsv ../../links.tsv
	mv multi.refine.fas ../../multi.fas
	mv pair.refine.fas ../../pair.fas
	cd ../..
	rm -rf Processing Pairwise Results Chr*.fa chr* *.sh
' ::: TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_rmmasked
cd ..

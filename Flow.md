```shell
mkdir "reference" && cd "reference" || exit
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
wget https://www.arabidopsis.org/download_files/Genes/Col-CEN%20genome%20assembly%20release/ColCEN.fasta

pigz -dc TAIR10_chr_all.fas.gz | faops split-name stdin .
rm ChrC.fa ChrM.fa
cat Chr*.fa >TAIR10_unmasked.fa
rm Chr*.fa
pigz -dc Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz | sed 's/>/>Chr/g' | faops split-name stdin .
rm ChrPt.fa ChrMt.fa
cat Chr*.fa >TAIR10_E58masked.fa
rm Chr*.fa
faops split-name ColCEN.fasta .
rm ChrC.fa ChrM.fa
cat Chr*.fa >ColCEN_unmasked.fa
rm Chr*.fa ColCEN.fasta
rm TAIR10_chr_all.fas.gz Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
```

```shell
mkdir "rm" && cd "rm" || exit
faops split-name ../TAIR10_unmasked.fa .
```

```shell
for i in 1 2 3 4 5; do
	RepeatMasker -spec arabidopsis -pa 16 -s -xsmall -e ncbi -dir . Chr${i}.fa
	RM2Bed.py -d . Chr${i}.fa.out
done
cat Chr*.fa_rm.bed >TAIR10.rm.fofn
cat TAIR10.rm.fofn | bedtools sort -i - >TAIR10_repeatmasker.out.bed
```

```shell
for i in 1 2 3 4 5; do
	trf Chr${i}.fa 2 7 7 80 10 50 15 -l 25 -h -ngs >Chr${i}.fa.dat
done
```

```python
import sys
import pandas as pd

input_dats = ["Chr1.fa.dat", "Chr2.fa.dat", "Chr3.fa.dat", "Chr4.fa.dat", "Chr5.fa.dat"]
output_bed = "TAIR10_trf.bed"
header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()

trf = []
for datf in input_dats:
    chrom = None
    sys.stderr.write("\r" + datf)
    with open(datf, 'r') as dat:
        for line in dat:
            splitline = line.split()
            if line.startswith("Sequence:"):
                chrom = int(line.split()[1].strip())
            elif line.startswith("@"):
                chrom = splitline[0][1:].strip()
            else:
                try:
                    int(splitline[0])
                except ValueError:
                    continue
                trf.append([chrom] + splitline[0: (len(header) - 1)])

trf_df = pd.DataFrame(trf, columns=header)
trf_df["start"] = trf_df["start"].astype(int)
trf_df.sort_values(by=["#chr", "start"], inplace=True)
trf_df.to_csv(output_bed, sep="\t", index=False)
```

```shell
cat TAIR10_trf.bed TAIR10_repeatmasker.out.bed \
	| cut -f 1-3 >>TAIR10.tmp.msk.bed
cut -f 1-3 TAIR10.tmp.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| awk '$3-$2 > 2000 {print $0}' \
	| bedtools merge -d 100 -i - >TAIR10.tmp2.msk.bed
cut -f 1-3 TAIR10.tmp.msk.bed TAIR10.tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../TAIR10_unmasked.fa >../TAIR10_rmmasked.fa
rm Chr* TAIR10_repeatmasker.out.bed TAIR10_trf.bed TAIR10.tmp.msk.bed TAIR10.tmp2.msk.bed TAIR10.rm.fofn
```

```shell
faops split-name ../ColCEN_unmasked.fa .
```

```shell
for i in 1 2 3 4 5; do
	RepeatMasker -spec arabidopsis -pa 16 -s -xsmall -e ncbi -dir . Chr${i}.fa
	RM2Bed.py -d . Chr${i}.fa.out
done
cat Chr*.fa_rm.bed >ColCEN.rm.fofn
cat ColCEN.rm.fofn | bedtools sort -i - >ColCEN_repeatmasker.out.bed
```

```shell
for i in 1 2 3 4 5; do
	trf Chr${i}.fa 2 7 7 80 10 50 15 -l 25 -h -ngs >Chr${i}.fa.dat
done
```

```python
import sys
import pandas as pd

input_dats = ["Chr1.fa.dat", "Chr2.fa.dat", "Chr3.fa.dat", "Chr4.fa.dat", "Chr5.fa.dat"]
output_bed = "ColCEN_trf.bed"
header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()

trf = []
for datf in input_dats:
    chrom = None
    sys.stderr.write("\r" + datf)
    with open(datf, 'r') as dat:
        for line in dat:
            splitline = line.split()
            if line.startswith("Sequence:"):
                chrom = int(line.split()[1].strip())
            elif line.startswith("@"):
                chrom = splitline[0][1:].strip()
            else:
                try:
                    int(splitline[0])
                except ValueError:
                    continue
                trf.append([chrom] + splitline[0: (len(header) - 1)])

trf_df = pd.DataFrame(trf, columns=header)
trf_df["start"] = trf_df["start"].astype(int)
trf_df.sort_values(by=["#chr", "start"], inplace=True)
trf_df.to_csv(output_bed, sep="\t", index=False)
```

```shell
cat ColCEN_trf.bed ColCEN_repeatmasker.out.bed \
	| cut -f 1-3 >>ColCEN.tmp.msk.bed
cut -f 1-3 ColCEN.tmp.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| awk '$3-$2 > 2000 {print $0}' \
	| bedtools merge -d 100 -i - >ColCEN.tmp2.msk.bed
cut -f 1-3 ColCEN.tmp.msk.bed ColCEN.tmp2.msk.bed \
	| bedtools sort -i - \
	| bedtools merge -i - \
	| seqtk seq -l 50 -M /dev/stdin ../ColCEN_unmasked.fa >../ColCEN_rmmasked.fa
rm Chr* ColCEN_repeatmasker.out.bed ColCEN_trf.bed ColCEN.tmp.msk.bed ColCEN.tmp2.msk.bed ColCEN.rm.fofn
```

```shell
cd .. && rm -rf "rm"
mkdir "smcover" && cd "smcover" || exit
egaz prepseq ../TAIR10_unmasked.fa -o .
faops masked ../TAIR10_unmasked.fa | spanr cover stdin | spanr stat --all chr.sizes stdin
faops masked ../TAIR10_rmmasked.fa | spanr cover stdin | spanr stat --all chr.sizes stdin
faops masked ../TAIR10_E58masked.fa | spanr cover stdin | spanr stat --all chr.sizes stdin
rm chr.* Chr.*
egaz prepseq ../ColCEN_unmasked.fa -o .
faops masked ../ColCEN_unmasked.fa | spanr cover stdin | spanr stat --all chr.sizes stdin
faops masked ../ColCEN_rmmasked.fa | spanr cover stdin | spanr stat --all chr.sizes stdin
cd .. && rm -rf "smcover"
cd ..
```

|                  | chrLength |     size | coverage |
|:-----------------|----------:|---------:|---------:|
| TAIR10_unmasked  | 119146348 |   186207 |   0.0016 |
| TAIR10_rmmasked  | 119146348 | 21324693 |   0.1790 |
| TAIR10_E58masked | 119146348 | 38123558 |   0.3200 |
| ColCEN_unmasked  | 131559676 |        0 |   0.0000 |
| ColCEN_rmmasked  | 131559676 | 33290733 |   0.2530 |

```shell
mkdir "biser" && cd "biser" || exit
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	biser -t 20 -o ${PREFIX}.out ../reference/${PREFIX}.fa
	awk '{print $1"("$9"):"$2"-"$3"\t"$4"("$10"):"$5"-"$6}' ${PREFIX}.out \
		| linkr sort stdin \
		| linkr clean stdin -o ${PREFIX}.links.sort.clean.tsv
done
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	linkr connect ${PREFIX}.links.sort.clean.tsv -r 0.05 \
		| linkr filter stdin -r 0.05 -o ${PREFIX}.links.filter.tsv
done
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	rgr merge ${PREFIX}.links.sort.clean.tsv -c 0.95 -o ${PREFIX}.links.merge.tsv
	linkr clean ${PREFIX}.links.sort.clean.tsv -r ${PREFIX}.links.merge.tsv --bundle 500 -o ${PREFIX}.links.clean.tsv
	linkr connect ${PREFIX}.links.clean.tsv -r 0.05 \
		| linkr filter stdin -r 0.05 -o ${PREFIX}.links.filter0.tsv
done
cd ..
```

```shell
echo "sample,key,count" >links.count.csv
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	for n in 2 3 4-50; do
		linkr filter ${PREFIX}.links.filter.tsv -n ${n} -o stdout >${PREFIX}.links.copy${n}.tsv
		wc -l ${PREFIX}.links.copy${n}.tsv \
			| perl -nl -e "
				@fields = grep { /\S+/ } split /\s+/;
				next unless @fields == 2;
				next unless \$fields[1] =~ /(\w+).links\.([\w-]+)\.tsv/;
				printf qq{%s,%s,%s\n}, \$1, \$2, \$fields[0];
        " \
				>>links.count.csv
	done
done
echo "sample,key,count" >links.count0.csv
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_unmasked ColCEN_rmmasked; do
	for n in 2 3 4-50; do
		linkr filter ${PREFIX}.links.filter0.tsv -n ${n} -o stdout >${PREFIX}.links.copy${n}.tsv
		wc -l ${PREFIX}.links.copy${n}.tsv \
			| perl -nl -e "
				@fields = grep { /\S+/ } split /\s+/;
				next unless @fields == 2;
				next unless \$fields[1] =~ /(\w+).links\.([\w-]+)\.tsv/;
				printf qq{%s,%s,%s\n}, \$1, \$2, \$fields[0];
        " \
				>>links.count0.csv
	done
done
```

```shell
rm *.links.merge.tsv *.links.sort.clean.tsv
cd ..
```

```shell
mkdir "lastz" && cd "lastz" || exit
echo 'strain,strain_id,species,species_id,genus,genus_id,family,family_id,order,order_id
Atha,3702,"Arabidopsis thaliana",3702,Arabidopsis,3701,Brassicaceae,3700,Brassicales,3699' >ensembl_taxon.csv
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_rmmasked; do
	mkdir ${PREFIX}
	egaz prepseq ../reference/${PREFIX}.fa -o ${PREFIX}/
	egaz template ${PREFIX}/ --self -o ${PREFIX}/ --taxon ./ensembl_taxon.csv --circos --parallel 16 -v
	cd ${PREFIX} || exit
	mkdir -p Pairwise
	egaz lastz \
		--isself --set set01 -C 0 \
		--parallel 16 --verbose \
		. . \
		-o Pairwise/AthavsSelf
	egaz lpcnam \
		--parallel 16 --verbose \
		. . \
		Pairwise/AthavsSelf
	cd ..
done
```

```shell
for PREFIX in TAIR10_unmasked TAIR10_rmmasked TAIR10_E58masked ColCEN_rmmasked; do
	cd ${PREFIX} || exit
	mkdir -p Results/Atha
	mkdir -p Processing/Atha
	ln -s "$(pwd)"/chr.fasta Processing/Atha/genome.fa
	cp -f "$(pwd)"/chr.sizes Processing/Atha/chr.sizes

	cd Processing/Atha || exit
	fasops axt2fas \
		../../Pairwise/AthavsSelf/axtNet/*.axt.gz \
		-l 1000 -s chr.sizes -o stdout >axt.fas
	fasops separate axt.fas -o . --nodash -s .sep.fasta
	egaz exactmatch target.sep.fasta genome.fa --length 500 --discard 50 -o replace.target.tsv
	fasops replace axt.fas replace.target.tsv -o axt.target.fas
	egaz exactmatch query.sep.fasta genome.fa --length 500 --discard 50 -o replace.query.tsv
	fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas
	fasops covers axt.correct.fas -o axt.correct.yml
	spanr split axt.correct.yml -s .temp.yml -o .
	spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
	spanr stat chr.sizes axt.union.yml -o union.csv
	fasops links axt.correct.fas -o stdout | perl -nl -e "s/(target|query)\.//g; print;" >links.lastz.tsv
	fasops separate axt.correct.fas --nodash --rc -o stdout \
		| perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" \
		| faops filter -u stdin stdout \
		| faops filter -n 250 stdin stdout \
			>axt.gl.fasta
	cd ../../..
done

```

```shell
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
```

| Name                                 | chrLength | defined size | coverage | copy2 count | BISER | BISER merged | ASGART |
|:-------------------------------------|:----------|:-------------|:---------|:------------|:------|:-------------|:-------|
| RepeatMasked with CONS-Dfam_withRBRM | 119146348 | 13098118     | 0.1096   | 3460        | 6314  | 2673         | 670    |


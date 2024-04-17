mkdir "lastz" && cd "lastz" || exit
echo 'strain,strain_id,species,species_id,genus,genus_id,family,family_id,order,order_id
Atha,3702,"Arabidopsis thaliana",3702,Arabidopsis,3701,Brassicaceae,3700,Brassicales,3699' >ensembl_taxon.csv
egaz prepseq ../TAIR10_masked.fa -o ./
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

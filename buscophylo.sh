cd /scratch//phylo
for i in `ls *.fasta |sed 's/.fasta//g'`;do busco -i $i".fna" -m genome -c 40 -l poales_odb10 --out $i; done
ls > namelist
#nano remove wrong names
for i in `cat namelist`; do ls /scratch/busco_sequences/single_copy_busco_sequences | grep 'fna' >$i.sbco.txt; done
sed -i 's/.fna//g' *.sbco.txt

module load r-env
start-r
Typha_latifolia<-scan("Typha_latifolia.sbco.txt", character(), quote = "")
Streptochaeta_angustifolia<-scan("Streptochaeta_angustifolia.sbco.txt",character(), quote = "")
Sorghum_bicolor<-scan("Sorghum_bicolor.sbco.txt",character(), quote = "")
Setaria_viridis<-scan("Setaria_viridis.sbco.txt",character(), quote = "")
Puya_raimondii<-scan("Puya_raimondii.sbco.txt",character(), quote = "")
Y7subA<-scan("Y7toref.subgenomeA.sbco.txt",character(), quote = "")
Y7subB<-scan("Y7toref.subgenomeB.sbco.txt",character(), quote = "")
Y21subA<-scan("Y21toref.subgenomeA.sbco.txt", character(), quote = "")
Y21subB<-scan("Y21toref.subgenomeB.sbco.txt", character(), quote = "")
Y17subA<-scan("Y17toref.subgenomeA.sbco.txt",character(), quote = "")
Y17subB<-scan("Y17toref.subgenomeB.sbco.txt",character(), quote = "")
EUsubA<-scan("PaEUtoref.subgenomeA.sbco.txt",character(), quote = "")
EUsubB<-scan("PaEUtoref.subgenomeB.sbco.txt",character(), quote = "")
Pharus_latifolius<-scan("Pharus_latifolius.sbco.txt",character(), quote = "")
Oryza_sativa<-scan("Oryza_sativa_japonica.sbco.txt",character(), quote = "")
Oropetium_thomaeum<-scan("Oropetium_thomaeum.sbco.txt",character(), quote = "")
Olyra_latifolia<-scan("Olyra_latifolia.sbco.txt",character(), quote = "")
#Mlu<-scan("Mlu_HiC_genome.sbco.txt",character(), quote = "")
Juncus_inflexus<-scan("Juncus_inflexus.sbco.txt",character(), quote = "")
Elaeis_guineensis<-scan("Elaeis_guineensis.sbco.txt",character(), quote = "")
Carex_cristatella<-scan("Carex_cristatella.sbco.txt",character(), quote = "")
Brachypodium_distachyon<-scan("Brachypodium_distachyon.sbco.txt",character(), quote = "")
Ananas_comosus<-scan("Ananas_comosus_F153.sbco.txt",character(), quote = "")
Aegilops_longissima<-scan("Aegilops_longissima.sbco.txt",character(), quote = "")
#Miscanthus_sinensis<-scan("Miscanthus_sinensis.sbco.txt",character(), quote = ""), this has too few sco compare to others, so remove
Cynodon_transvaalensis<-scan("Cynodon_transvaalensis.sbco.txt",character(), quote = "")
Echinochloa_haploclada<-scan("Echinochloa_haploclada.sbco.txt",character(), quote = "")
Pa_CNsubgenomeA<-scan("CN.subgenomeA.sbco.txt",character(), quote = "")
Pa_CNsubgenomeB<-scan("CN.subgenomeB.sbco.txt",character(), quote = "")
Amborellatri<-scan("Amborellatri.sbco.txt",character(), quote = "")
d<-Reduce(intersect, list(Y7subA,Y7subB,Y21subA,Y21subB,Y17subA,Y17subB,EUsubA, EUsubB,Pa_CNsubgenomeA,Pa_CNsubgenomeB,Typha_latifolia,Streptochaeta_angustifolia,Sorghum_bicolor,Setaria_viridis,Puya_raimondii,Pharus_latifolius,Oryza_sativa,Oropetium_thomaeum,Olyra_latifolia,Juncus_inflexus,Elaeis_guineensis,Carex_cristatella,Amborellatri,Brachypodium_distachyon,Ananas_comosus,Aegilops_longissima,Cynodon_transvaalensis,Echinochloa_haploclada))
write.csv(d,"commonsbco_subgenome.txt", quote=F, row.names=FALSE)
#remove row name
sed -i '1d' commonsbco_subgenome.txt
grep -w -f commonsbco.txt 
#do it one by one
#cat for each busco of different species
for i in `cat commonsbco_subgenome.txt `; do for j in `cat namelist`; do cat $j/run_poales_odb10/busco_sequences/single_copy_busco_sequences/$i.fna | sed 's/^>.*$/>'"$j"'/g' >> $i.all.fna; done; done
#mafft align
for i in *.all.fna; do mafft $i > $i.aligned; done
module load raxml
for i in `ls *.aligned| sed 's/\.all.fna.aligned//g'`; do raxmlHPC -f a -s $i.all.fna.aligned -x 1245 -N 100 -m GTRGAMMA -n $i -T 8 -p 5697; done
#astral 5.7.3
java -jar  Astral/astral.5.7.8.jar 
for i in RAxML_bestTree.*; do cat $i >> astral.tree; done
ls $PWD/RAxML_bootstrap.* > astral.BS.txt
mv astral* ..
cd ..
java -jar  ./ASTRAL/Astral/astral.5.7.8.jar -i astral.tree -b astral.BS.txt -o species.bootstraped.astral.tre
#divergence time estimation using mcmctree,cyrus
cd /data2/mcmctree

for i in *.all.fna.aligned; do seqkit seq -w 0 $i > $i.oneline ; done
for i in *.all.fna.aligned.oneline; do awk 'FNR==2{print FILENAME,length}' $i >>OG.ID; done
sed -i 's/.all.fna.aligned.oneline//g' OG.ID
mv OG.ID OG.length
cat ./OG.length | while read line
do
        sample_id=$(echo $line |awk '{print $1}') #get OG.ID
        sample_a=$(echo $line |awk '{print $2}') #get sequence length
        cat ${sample_id}.all.fna.aligned.oneline |seqkit seq -w 0|sed -E ":a;N;s/\n/ /g;ta" |sed "s/ >/\n/g" |sed "s/>//g"|sed "s/ /  /g"|sed "1i\28  ${sample_a}" |sed '1i\ ' > ./${sample_id}.phy #convert the ${sample_id}.fa obtained from the last steps to phylip format，and add empty line at the initial lines so that each OG can be separated after being merged. DO remember to change the number of species, here is 28

done
cd ..
cat ./phylo/*phy >input.phy 
#use Amborella as outgroup, manual change setting files, and use mcmctree


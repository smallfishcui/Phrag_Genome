#cyrus
source activate ragtag
#ln -s /data2/gapclose/Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta .
ln -s /data2/genomephasing/PaEUcurated.fasta .
ln -s /data2/genome/refY17/Y17_curated.fasta .
ln -s /data2/genome/refY21/Y21_curated.fasta .
ln -s /data2/genome/refY7/Y7_curated.fasta .
ragtag.py scaffold Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta PaEUcurated.fasta -C -u -r
mv ragtag_output PaEUtoref
ragtag.py scaffold Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta Y21_curated.fasta -C -u -r
mv ragtag_output Y21toref
ragtag.py scaffold Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta Y17_curated.fasta -C -u -r
mv ragtag_output Y17toref
ragtag.py scaffold Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta Y7_curated.fasta -C -u -r
mv ragtag_output Y7toref

#lifoff gene nodels
cd /data2/ragtag
ln -s /data2/gapclose/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
conda activate liftoff
liftoff -g EVM.all.gff3 ./Y21toref.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y21toref/Y21.liftoff.gff
liftoff -g EVM.all.gff3 ./Y17toref.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y17.liftoff.gff
liftoff -g EVM.all.gff3 ./Y7toref.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y7.liftoff.gff
liftoff -g EVM.all.gff3 ./PaEUtoref.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Pa_EU.liftoff.gff

cd /data2/liftoff
ln -s /data2/gapclose/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
ln -s /data2/ragtag/EVM.all.gff3 .
ln -s /data2/genomephasing/PaEUcurated.fasta .
ln -s /data2/genome/refY17/Y17_curated.fasta .
ln -s /data2/genome/refY21/Y21_curated.fasta .
ln -s /data2/genome/refY7/Y7_curated.fasta .
ln -s /data2/ragtag/PaEUtoref.fasta .
liftoff -g EVM.all.gff3 ./Y21_curated.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y21_curated.liftoff.gff
liftoff -g EVM.all.gff3 ./Y17_curated.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y17_curated.liftoff.gff
liftoff -g EVM.all.gff3 ./Y7_curated.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Y7_curated.liftoff.gff
liftoff -g EVM.all.gff3 ./PaEUcurated.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./PaEUcurated.liftoff.gff
liftoff -g PaEU_purged.gff ./PaEUtoref.fasta PaEUcurated.fasta -o ./PaEUtoref.liftoff.original.gff
/opt/EVidenceModeler/EvmUtils/gff3_file_to_proteins.pl ./Y21/Y21_curated.liftoff.gff Y21_curated.fasta > Y21.faa

#gemoma does better than liftoff,use gemoma instead to get the annotation for each original assemblies
ln -s /scratch/genome/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
module load openjdk/11.0.12_7
export PATH=/scratch/gemoma/mmseqs/bin:$PATH
java -jar -Xmx200000m /scratch/gemoma/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=32 outdir=./PaEU_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=ragtag.scaffold.fasta a=EVM.all.gff3 g=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta &> Pa_Osat_ann.log
perl ./GeMoMa_gff_to_gff3.pl ./final_annotation.gff > PaEU_gemoma_evm_mod.gff3 # gene number 42006, much less than its de novo annotation. use its own annotation
#use gemoma also for Y7, Y21, and Y17
cd /scratch/ragtag
module load openjdk/11.0.12_7
export PATH=/scratch/gemoma/mmseqs/bin:$PATH
java -jar -Xmx200000m /scratch/gemoma/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=32 outdir=./Y17_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Y17_curated.fasta a=EVM.all.gff3 g=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta &> Y17_ref_ann.log#gene number 40986, more than liftoff 40688
java -jar -Xmx200000m /scratch/gemoma/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=32 outdir=./Y21_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Y21_curated.fasta a=EVM.all.gff3 g=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta &> Y21_ref_ann.log#40232, less than liftoff 40969
java -jar -Xmx200000m /scratch/gemoma/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=32 outdir=./Y7_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Y7_curated.fasta a=EVM.all.gff3 g=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta &> Y7_ref_ann.log#gene number 38903,less than liftoff 40124

#LAI
conda activate genometools
source activate ltr_retriever
cd /data2/LAI
mkdir PaEUtoref
cd PaEUtoref
ln -s /data2/ragtag/PaEUtoref.fasta .
grep 'PaChr' PaEUtoref.fasta | sed 's/>//g' > chrom.list
samtools faidx PaEUtoref.fasta -r chrom.list > PaEUtoref.chronly.fasta
sed 's/RagTag/PaEU/g' PaEUtoref.chronly.fasta > PaEUtoref.renamed.fasta
sed -i 's/HiC_scaffold/Chr/g' PaEUtoref.renamed.fasta
gt suffixerator -db PaEUtoref.renamed.fasta -indexname PaEUtoref.renamed.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index PaEUtoref.renamed.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > PaEUtoref.renamed.fasta.harvest.scn
/data2/LAI/EDTA/bin/LTR_FINDER_parallel/LTR_FINDER_parallel -seq PaEUtoref.renamed.fasta -threads 10 -harvest_out -size 1000000 -time 300
cat PaEUtoref.renamed.fasta.harvest.scn PaEUtoref.renamed.fasta.finder.combine.scn > PaEUtoref.renamed.fasta.rawLTR.scn
/data2/LAI/LTR_retriever/LTR_retriever -genome PaEUtoref.renamed.fasta -inharvest PaEUtoref.renamed.fasta.rawLTR.scn -threads 10 
cd ..

mkdir Y21toref
cd Y21toref
ln -s /data2/ragtag/Y21toref.fasta .
grep 'PaChr' Y21toref.fasta | sed 's/>//g' > chrom.list
samtools faidx Y21toref.fasta -r chrom.list > Y21toref.chronly.fasta
sed 's/RagTag/Y21/g' Y21toref.chronly.fasta > Y21toref.renamed.fasta
sed -i 's/HiC_scaffold/Chr/g' Y21toref.renamed.fasta
gt suffixerator -db Y21toref.renamed.fasta -indexname Y21toref.renamed.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Y21toref.renamed.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > Y21toref.renamed.fasta.harvest.scn
/data2/LAI/EDTA/bin/LTR_FINDER_parallel/LTR_FINDER_parallel -seq Y21toref.renamed.fasta -threads 10 -harvest_out -size 1000000 -time 300
cat Y21toref.renamed.fasta.harvest.scn Y21toref.renamed.fasta.finder.combine.scn > Y21toref.renamed.fasta.rawLTR.scn
/data2/LAI/LTR_retriever/LTR_retriever -genome Y21toref.renamed.fasta -inharvest Y21toref.renamed.fasta.rawLTR.scn -threads 10 
cd ..

mkdir Y17toref
cd Y17toref
ln -s /data2/ragtag/Y17toref.fasta .
grep 'PaChr' Y17toref.fasta | sed 's/>//g' > chrom.list
samtools faidx Y17toref.fasta -r chrom.list > Y17toref.chronly.fasta
sed 's/RagTag/Y17/g' Y17toref.chronly.fasta > Y17toref.renamed.fasta
sed -i 's/HiC_scaffold/Chr/g' Y17toref.renamed.fasta
gt suffixerator -db Y17toref.renamed.fasta -indexname Y17toref.renamed.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Y17toref.renamed.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > Y17toref.renamed.fasta.harvest.scn
/data2/LAI/EDTA/bin/LTR_FINDER_parallel/LTR_FINDER_parallel -seq Y17toref.renamed.fasta -threads 10 -harvest_out -size 1000000 -time 300
cat Y17toref.renamed.fasta.harvest.scn Y17toref.renamed.fasta.finder.combine.scn > Y17toref.renamed.fasta.rawLTR.scn
/data2/LAI/LTR_retriever/LTR_retriever -genome Y17toref.renamed.fasta -inharvest Y17toref.renamed.fasta.rawLTR.scn -threads 10 
cd ..

mkdir Y7toref
cd Y7toref
ln -s /data2/ragtag/Y7toref.fasta .
grep 'PaChr' Y7toref.fasta | sed 's/>//g' > chrom.list
samtools faidx Y7toref.fasta -r chrom.list > Y7toref.chronly.fasta
sed 's/RagTag/Y7/g' Y7toref.chronly.fasta > Y7toref.renamed.fasta
sed -i 's/HiC_scaffold/Chr/g' Y7toref.renamed.fasta
gt suffixerator -db Y7toref.renamed.fasta -indexname Y7toref.renamed.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Y7toref.renamed.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > Y7toref.renamed.fasta.harvest.scn
/data2/LAI/EDTA/bin/LTR_FINDER_parallel/LTR_FINDER_parallel -seq Y7toref.renamed.fasta -threads 10 -harvest_out -size 1000000 -time 300
cat Y7toref.renamed.fasta.harvest.scn Y7toref.renamed.fasta.finder.combine.scn > Y7toref.renamed.fasta.rawLTR.scn
/data2/LAI/LTR_retriever/LTR_retriever -genome Y7toref.renamed.fasta -inharvest Y7toref.renamed.fasta.rawLTR.scn -threads 10 
cd ..

mkdir ref
cd ref
ln -s /data2/gapclose/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
grep 'PaChr' Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta | sed 's/>//g' > chrom.list
samtools faidx Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -r chrom.list > Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta
gt suffixerator -db Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta -indexname Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.harvest.scn
/data2/LAI/EDTA/bin/LTR_FINDER_parallel/LTR_FINDER_parallel -seq Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta -threads 10 -harvest_out -size 1000000 -time 300
cat Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.harvest.scn Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.finder.combine.scn > Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.rawLTR.scn
/data2/LAI/LTR_retriever/LTR_retriever -genome Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta -inharvest Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.rawLTR.scn -threads 10 
cd ..

#create a manhattan plot for LAI
cd /data2/LAI/ref
awk '{print $1"\t"$3"\t"$7}' Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.out.LAI | grep 'PaChr' > Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.out.chronly.mahattanplot.LAI
#own computer, R
setwd("/Users/Downloads")
library(data.table)
library(qqman)
a<-read.delim("Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.renamed.fasta.out.chronly.mahattanplot.LAI", header=T)
setDT(a, keep.rownames = TRUE)[]
manhattan(a, chr="Chr", bp="to", snp="rn", p="LAI", logp = FALSE, suggestiveline = FALSE, genomewideline = 10)


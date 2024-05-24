cd /data2/orthofinder
/opt/EVidenceModeler/EvmUtils/gff3_file_to_proteins.pl PaEUtoref.liftoff.original.gff PaEUtoref.fasta > PaEUtoref.faa
orthofinder -f ./ -a 40 
cd /scratch/orthofinder
start-r
.libPaths(c("/projappl/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
if(!requireNamespace('BiocManager', quietly = TRUE))
install.packages('BiocManager')
BiocManager::install("cogeqc")
library(ggtree)
library(cogeqc)
library(treeio)
stats_dir <-"./Comparative_Genomics_Statistics"
ortho_stats <- read_orthofinder_stats(stats_dir)
tree<-read.newick("SpeciesTree_rooted.txt")
p1<-plot_species_tree(tree, xlim = c(-1, 1), stats_list = ortho_stats)
p2<-plot_duplications(ortho_stats)
p3<-plot_genes_in_ogs(ortho_stats)
p4<-plot_species_specific_ogs(ortho_stats)
pdf("duplications.pdf")
p2
dev.off()
pdf("genes_in_ogs.pdf")
p3
dev.off()
pdf("species_specific_ogs.pdf")
p4
dev.off()

#invasive lineage specific orthogroups
cd /data2/OrthoFinder/Results_May11/Orthogroups
awk -F'\t' '$2==0 && $4==0 && $5==0 && $6==0 { print $0 } ' Orthogroups.GeneCount.tsv | awk '{ print $1 }' > PaEU_specific.orthogroup
cd /data2/OrthoFinder/Results_May11/Orthogroup_Sequences
find . -name "*.fa" | grep -f /data2/OrthoFinder/Results_May11/Orthogroups/PaEU_specific.orthogroup | xargs cat > Pau.specific.Orthogroups.fa
grep '>' Pau.specific.Orthogroups.fa | sed 's/>//g' > Pau.specific.Orthogroups.list

#test gene overlap between paEu and tandem duplications
grep '>' PaEUtoref.faa | awk '{ print $1 }' | sed 's/>//g' > all.gene
comm -12 Pau.specific.Orthogroups.sorted.list invasitoinvasi.tandems.sorted.genelist > comm.genelist
start-r
.libPaths(c("/projappl/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(GeneOverlap)
a<-read.table("Pau.specific.Orthogroups.sorted.list",header=F)
b<-read.table("invasitoinvasi.tandems.sorted.genelist",header=F)
allgene<-read.table("all.gene", header=F)
go.obj <- newGeneOverlap(a$V1,b$V1,genome.size=length(allgene$V1))
go.obj <- testGeneOverlap(go.obj)
go.obj

#This script provided the pipeline for assembling a plant genome using pacbio HiFi reads and scaffolding with HiC contacts. 
#It is NOT an automated pipeline! Folder structure is not well illustrated because several steps were performed on different servers.
#It very important that you know about the meaning of each software. Manual correction needs to be incorporated in the steps, so check each step while using this pipeline.
#We first assemble the chloroplast and mitochondria genome assembly, and then remove the reads that are mapped to these organelles.

#assemble chloroplast genome assembly using the HiFi reads. Note, mitochondria genome can also be obtained here, however, it is more difficult to get the circular sequences than flye.
cat P05TYD22610974-1_r64400e_20220516_080759_1_C02.hifi_reads.fasta P05TYD22610974-1_r64400e_20220519_071347_1_E02.hifi_reads.fasta > combined.fasta
PMAT autoMito -i combined.fasta -o ./Pamt -st hifi -g 920m -m -cpu 40 -fc 0.2 
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -c ./Pamt/assembly_result/PMATContigGraph.txt -a ./Pamt/assembly_result/PMATAllContigs.fna -gs 920m -rs ./Pamt/subsample/assembly_seq.cut20K.fasta -tp pt -o ./chloroplast_gfa -s 343 345 905 513 1344
cd ./chloroplast_gfa/gfa_result
awk '/^S/{print ">"$2"\n"$3}' PMAT_pt_master.gfa > PMAT_pt_master.fa #manually curate using Bandage
#awk '/^S/{print ">"$2"\n"$3}' PMAT_mt_master.gfa > PMAT_mt_master.fa

#assemble mitochondria genome
source activate flye
flye --meta --pacbio-hifi ./combined.fasta -o ./ --threads 40  
get_organelle_from_assembly.py -F embplant_mt -g assembly_graph.gfa -t 48 -o ./mito --min-depth 100 --max-depth 3000 --overwrite 
awk '/^S/{print ">"$2"\n"$3}' embplant_mt.complete.graph1.selected_graph.gfa > embplant_mt.complete.graph1.selected_graph.fasta

#Take the reads that were unmapped to the organelles and use them to assemble the nuclear genome
cd ./HiFireadsnoorganelle
ln -s ../PAMT/combined.fasta .
#get the chloroplast and mitochondria assembly and place them in the same folder
cat Chloroplast_path_sequence.fasta embplant_mt.complete.graph1.1.path_sequence.fasta > organelle.fasta
minimap2 -t 40 --secondary=no -ax map-hifi ./organelle.fasta combined.fasta -a | samtools view -f 4 | samtools fastq -@ 40 > reads.HiFinoorganelle.fastq2

#nevigate to a new folder
cd ../Hifiasm
ln -s ../HiFireadsnoorganelle/reads.HiFinoorganelle.fastq .
./hifiasm/hifiasm -o Pa_noorganelle.asm -t40 reads.HiFinoorganelle.fastq
awk '/^S/{print ">"$2;print $3}' Pa_noorganelle.asm.bp.p_ctg.gfa > Pa_noorganelle.asm.bp.p_ctg.fa

#use allhic pipeline, first build synteny with reference species Oryza sativa
source activate gemoma
cd ../gemoma
ln -s ../hifiasm/Pa_noorganelle.asm.bp.p_ctg.fa .
GeMoMa GeMoMaPipeline threads=24 outdir=./phraghifi_noorganelpredic GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Pa_noorganelle.asm.bp.p_ctg.fa a=GCF_001433935.1_IRGSP-1.0_genomic.gff g=GCF_001433935.1_IRGSP-1.0_genomic.fna -Xmx6000m &> phraghifipredic.log

#start allhic pipeline, LW_S0_L000_R1_000.fastq.gz and LW_S0_L000_R2_000.fastq.gz are the hic fastq files
cd ../allhic
bwa index -a bwtsw Pa_noorganelle.asm.bp.p_ctg.fa
samtools faidx Pa_noorganelle.asm.bp.p_ctg.fa
/opt/bwa/bwa aln -t 20 Pa_noorganelle.asm.bp.p_ctg.fa LW_S0_L000_R1_000.fastq.gz > sample_R1.sai  
/opt/bwa/bwa aln -t 40 Pa_noorganelle.asm.bp.p_ctg.fa LW_S0_L000_R2_000.fastq.gz > sample_R2.sai 
/opt/bwa/bwa sampe Pa_noorganelle.asm.bp.p_ctg.fa sample_R1.sai sample_R2.sai LW_S0_L000_R1_000.fastq.gz LW_S0_L000_R2_000.fastq.gz > sample.bwa_aln.sam
PreprocessSAMs.pl sample.bwa_aln.sam Pa_noorganelle.asm.bp.p_ctg.fa MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -bt Pa_noorganelle.asm.bp.p_ctg.fa.fai sample.clean.sam > sample.clean.bam 

#used Allele.ctg.table to prune 1) signals that link alleles and 2) weak signals from BAM files
ln -s ../gemoma/phraghifi_noorganelpredic/final_annotation.gff .
cat final_annotation.gff | grep 'mRNA'  > gemoma.gene.gff
bedtools getfasta -fi Pa_noorganelle.asm.bp.p_ctg.fa -bed  gemoma.gene.gff -fo Pa_noorganelle.asm.bp.p_ctg.gene.fa
#this line of code has been modified from the previous versions of coffee, birch, check carefully about the content of each column
awk '{ print $1":"$4-1"-"$5"\t"$9 }' gemoma.gene.gff | awk -F";" '{ print $1"\t"$2}'| awk -F"\t" '{ print $1"\t"$2}' | sed 's/ID=//g' > gemoma.map.txt
awk 'NR==FNR{a[">"$1]=">"$2; next} /^>/{$0=a[$0]} 1' gemoma.map.txt Pa_noorganelle.asm.bp.p_ctg.gene.fa  > Pa_noorganelle.asm.bp.p_ctg.gene.rename.fa
awk -F "\t" '$3=="mRNA"{ print $0 }' GCF_001433935.1_IRGSP-1.0_genomic.gff > Oryza_sativa.gene.gff

#remember to delete old index file for the Osativa genome, otherwise it doesnt work
bedtools getfasta -fi GCF_001433935.1_IRGSP-1.0_genomic.fna -bed Oryza_sativa.gene.gff -fo Oryza_sativa.gene.fasta

#modified from previous scripts
cat Oryza_sativa.gene.gff | awk '{ print $1":"$4-1"-"$5"\t"$9 }' | awk -F";" '{ print $1"\t"$2 }' | awk -F"\t" '{ print $1"\t"$3 }'| sed 's/Parent=//g' > Oryza_sativa.map.txt
awk 'NR==FNR{map[">"$1]=$2; next} /^>/{$0=">"map[$0]} 1' Oryza_sativa.map.txt Oryza_sativa.gene.fasta > Oryza_sativa.rename.fasta
makeblastdb -in Oryza_sativa.rename.fasta -dbtype nucl
blastn -query Pa_noorganelle.asm.bp.p_ctg.gene.rename.fa -db Oryza_sativa.rename.fasta -out hifi_vs_Osativa.blast.out -evalue 0.001 -outfmt 6 -num_threads 40 -num_alignments 1
blastn_parse.pl -i hifi_vs_Osativa.blast.out -o Ehifi_vs_final.blast.out -q Pa_noorganelle.asm.bp.p_ctg.gene.rename.fa -b 1 -c 0.6 -d 0.8 
sed -i 's/Parent/Name/g' Oryza_sativa.gene.gff
sed -i 's/mRNA/gene/g' Oryza_sativa.gene.gff
sed -i 's/ID/Name/g' gemoma.gene.gff
sed -i 's/mRNA/gene/g' gemoma.gene.gff
classify.pl -i Ehifi_vs_final.blast.out -p 2 -r Oryza_sativa.gene.gff -g gemoma.gene.gff
partition.pl -g Allele.gene.table -r Pa_noorganelle.asm.bp.p_ctg.fa -b sample.clean.bam -d wrk_dir

cd ./wrk_dir
for i in NC_*; do 
cd $i;
ALLHiC_prune -i ../../Allele.ctg.table -b sample.clean.bam -r seq.fasta
ALLHiC_partition -r seq.fasta -b prunning.bam -e GATC -k 2
ALLHiC_rescue -r seq.fasta -b sample.clean.bam -c prunning.clusters.txt -i prunning.counts_GATC.txt
allhic extract sample.clean.bam seq.fasta --RE GATC
allhic optimize sample.clean.counts_GATC.txt sample.clean.clm
ALLHiC_build seq.fasta
cd ..;
done

#global rescue
cd wrk_dir
find * -maxdepth 0|while read line;do
agp=$line/groups.agp
grep -v '#' ${agp} | egrep 'sample.clean.counts_GATC' \
    |cut -f 1 |sort -u| while read chrID;do
    ctgNum=`grep -w ${chrID} ${agp} | grep 'W' |wc -l`;
    ctg=`grep -w ${chrID} ${agp} | grep -w  'W' | cut -f 6 | tr '\n' ' '`;
    echo -e "${line}\t${ctgNum}\t${ctg}" >> clusters.txt;
    done
done
cd ..

#global rescue
allhic extract sample.clean.bam Pa_noorganelle.asm.bp.p_ctg.fa --RE GATC
ALLHiC_rescue -r Pa_noorganelle.asm.bp.p_ctg.fa -b sample.clean.bam -c ./wrk_dir/clusters.txt -i sample.clean.counts_GATC.txt
for K in {1..14}; do allhic optimize group${K}.txt sample.clean.clm;done
ALLHiC_build Pa_noorganelle.asm.bp.p_ctg.fa
seqkit sort -l -r groups.asm.fasta > groups.asm.sorted.fasta
perl ./getFaLen.pl -i groups.asm.sorted.fasta -o len.txt
grep 'group' len.txt > chrn.list
/usr/bin/quast.py -s -e groups.asm.sorted.fasta
python ./ALLHiC/bin/ALLHiC_plot -b sample.clean.bam -a groups.agp -l chrn.list -s 500k

#convert allhic results to juicebox
samtools sort -n -@ 8 sample.bwa_aln.sam  -o sample.sort.bam
source activate matlock
matlock bam2 juicer sample.sort.bam out.links.txt
# sort out.links.txt
sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
#scripts from 3ddna github page
python agp2assembly.py groups.agp out.assembly
# 3d-dna to generate .hic
cd /data2/Cui/Phragmites/CNWGS/Genome_hifiasm_V4/allhic
#manually check hic file with juicebox
#convert assembly file to agp
python juicebox_assembly_converter.py -a Panoorganelle.review.firstround.assembly -f Pa_noorganelle.asm.bp.p_ctg.fa 

#manually adjust the assembly file based on telomere information and allhic results.
#check if the assembly works
python juicebox_assembly_converter.py -a Panoorganelle.review.firstround.revised.assembly -f Pa_noorganelle.asm.bp.p_ctg.fa 

#use juicer and 3ddna scaffold with HiC contacts, make the folder structure as the pipeline instructed, prepare the files
cd ../juicer
mkdir restriction_sites
python generate_site_positions.py DpnII Pa_noorganelle.asm.bp.p_ctg.fa ./references/Pa_noorganelle.asm.bp.p_ctg.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' Pa_noorganelle.asm.bp.p_ctg.fa_DpnII.txt > Pa_noorganelle.asm.bp.p_ctg.fa.chrom.sizes
export PATH="/opt/bwa/:$PATH"
./scripts/juicer.sh -g Pa_noorganelle.asm.bp.p_ctg.fa -z ./references/Pa_noorganelle.asm.bp.p_ctg.fa -D . -s DpnII -p ./restriction_sites/Pa_noorganelle.asm.bp.p_ctg.fa.chrom.sizes -y ./restriction_sites/Pa_noorganelle.asm.bp.p_ctg.fa_DpnII.txt -t 40 &> juicer.log
cd ..

#python agp2assembly.py Panoorganelle.review.firstround.agp Panoorganelle.review.firstround.revised.assembly
bash /opt/3d-dna/visualize/run-assembly-visualizer.sh -q 1 -p true Panoorganelle.review.firstround.revised.assembly  out.sorted.links.txt
awk -f ./wrap-fasta.awk Pa_noorganelle.asm.bp.p_ctg.fa > Pa_noorganelle.asm.bp.p_ctg.wrapped.fa
bash /opt/3d-dna/run-asm-pipeline-post-review.sh --sort-output -r Pa_noorganelle.asm.bp.p_ctg.wrapped.final.review.assembly Pa_noorganelle.asm.bp.p_ctg.wrapped.fa merged_nodups.txt  &> postreview.log









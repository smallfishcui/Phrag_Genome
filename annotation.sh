#mahti
singularity run dfam-tetools-latest.sif BuildDatabase -engine ncbi -name phrag Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta
singularity run dfam-tetools-latest.sif RepeatModeler -database phrag -engine ncbi -threads 128 
#make sure you softlink the classified file, otherwise you will not get a table of classified elements after the run.
ln -s ./RM_3026906.WedApr170825472024/consensi.fa.classified .
#softmaskinng
singularity run dfam-tetools-latest.sif  RepeatMasker -pa 32 -xsmall -lib consensi.fa.classified Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta

#gemoma, mahti
cd ./gemoma
java -jar -Xmx200000m GeMoMa-1.9.jar CLI GeMoMaPipeline threads=128 outdir=./Pa_annot_Osativ GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=./rice/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta a=./rice/GCF_001433935.1_IRGSP-1.0_genomic.gff g=./rice/GCF_001433935.1_IRGSP-1.0_genomic.fna r=MAPPED ERE.m=/scratch/project_2009273/RNAseq/merged_20.bam &> Pa_Osat_ann.log
srun java -jar -Xmx200000m GeMoMa-1.9.jar CLI GeMoMaPipeline threads=128 outdir=./Pa_annot_Phalli GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=./Phalli/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta a=./Phalli/Phallii_590_v3.2.gene.gff3 g=./Phalli/Phallii_590_v3.0.fa r=MAPPED ERE.m=/scratch/project_2009273/RNAseq/merged_20.bam &> Pa_Phalli_ann.log
srun java -jar -Xmx200000m GeMoMa-1.9.jar CLI GeMoMaPipeline threads=128 outdir=./Pa_annot_Sbicolor GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=./Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta a=./sorghum/Sbicolor_454_v3.1.1.gene.gff3 g=./sorghum/Sbicolor_454_v3.0.1.fa r=MAPPED ERE.m=/scratch/project_2009273/RNAseq/merged_20.bam &> Pa_sbicolor_ann.log
srun java -jar -Xmx200000m GeMoMa-1.9.jar CLI GeMoMaPipeline threads=128 outdir=./Pa_annot_Sviridis GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=./Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta a=./Sviridis/Sviridis_311_v1.1.gene.gff3 g=./Sviridis/Sviridis_311_v1.0.fa r=MAPPED ERE.m=/scratch/project_2009273/RNAseq/merged_20.bam &> Pa_Sviridis_ann.log
#GeMoMa_gff_to_gff3.pl downloaded from EVM folder
perl ./GeMoMa_gff_to_gff3.pl ./Pa_annot_Osativ/final_annotation.gff > Osativa_gemoma_evm_mod.gff3
perl ./GeMoMa_gff_to_gff3.pl ./Pa_annot_Phalli/final_annotation.gff > Phalli_gemoma_evm_mod.gff3
perl ./GeMoMa_gff_to_gff3.pl ./Pa_annot_Sbicolor/final_annotation.gff > Sbicolor_gemoma_evm_mod.gff3
perl ./GeMoMa_gff_to_gff3.pl ./Pa_annot_Sviridis/final_annotation.gff > Sviridis_gemoma_evm_mod.gff3

#annotation genemark-es, puhti, dont't use mahti
module load bioperl
cd ../annotation
./gmes_linux_64_4/gmes_petap.pl --ES --cores 40 --sequence Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta
perl /opt/EVidenceModeler/EvmUtils/misc/genemark_gtf2gff3.pl genemark.gtf > genemark_evm_mod.gff3 
mv genemark_evm_mod.gff3  genemark_es_evm_mod.gff3

#braker, puhti
module load biokit
cd /scratch/annotation/RNAseq
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p name.list)
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir ./ --genomeFastaFiles Pa_noorganelle.asm.bp.p_ctg.unwrapped.FINAL.sorted.gapclose.fasta.masked --genomeSAindexNbases 13
mkdir $name;
cd $name;
STAR --runMode alignReads --runThreadN 40 --genomeDir ../ --readFilesIn ../$name"_1.fastq.gz" ../$name"_2.fastq.gz" --outFilterMismatchNmax 15 --readFilesCommand zcat --outSAMstrandField intronMotif 
samtools sort -@ 40 Aligned.out.sam -o Aligned.out.bam.sort -O BAM
awk '{ print "./"$1"/Aligned.out.bam.sort" }' name.list > bam.list
samtools merge -o merged_20.bam -b bam.list -@40

#mahti
cd /scratch/braker
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
gunzip Viridiplantae.fa.gz
git clone https://github.com/Gaius-Augustus/Augustus.git
#use script
export PATH=/scratch/braker/Augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/scratch/braker/Augustus/config/
ln -s /scratch/repeat/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.masked .
ln -s /scratch/RNAseq/*.fastq.gz .
paste -sd, name.list > name.comma.list
singularity exec -B ${PWD}:${PWD} --env AUGUSTUS_CONFIG_PATH=/scratch/braker/Augustus/config/ braker3.sif braker.pl --species=Paustralis --genome=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.masked \
       --rnaseq_sets_ids=R20033279-P200805-SD-S136-1,R20033279-P200805-SD-S15-1,R20033279-P200805-SD-S150-1,R20033279-P200805-SD-S162-1,R20033279-P200805-SD-S174-1,R20033279-P200805-SD-S188-1,R20033279-P200805-SD-S188-2,R20033279-P200805-SD-S191-1,R20033279-P200805-SD-S207-1,R20033279-P200805-SD-S620-1,R20033279-P200805-SD-S68-1,R20033279-P200805-SD-lanzhou-1,R20033279-P200805-SD-yuncheng-1,SRR16648326,SRR16648327,SRR16648328,SRR16648329,SRR16648330,SRR16648331,SRR16648332,Unknown_BA854-01T0001,Unknown_BA854-01T0002,Unknown_BA854-01T0003 --rnaseq_sets_dirs=./ --threads=$OMP_NUM_THREADS --busco_lineage=poales_odb10
singularity exec -B ${PWD}:${PWD} --env AUGUSTUS_CONFIG_PATH=/scratch/braker/Augustus/config/ braker3.sif braker.pl --genome=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.masked \
      --prot_seq=./Viridiplantae.fa --threads=$OMP_NUM_THREADS --busco_lineage=poales_odb10
#etp
singularity exec -B /scratch/RNAseq/:/scratch/RNAseq/ -B ${PWD}:${PWD} --env AUGUSTUS_CONFIG_PATH=/scratch/braker/Augustus/config/ braker3.sif braker.pl --species=Paustralis2 --genome=Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.masked \
       --rnaseq_sets_ids=R20033279-P200805-SD-S136-1,R20033279-P200805-SD-S15-1,R20033279-P200805-SD-S150-1,R20033279-P200805-SD-S162-1,R20033279-P200805-SD-S174-1,R20033279-P200805-SD-S188-1,R20033279-P200805-SD-S188-2,R20033279-P200805-SD-S191-1,R20033279-P200805-SD-S207-1,R20033279-P200805-SD-S620-1,R20033279-P200805-SD-S68-1,R20033279-P200805-SD-lanzhou-1,R20033279-P200805-SD-yuncheng-1,SRR16648326,SRR16648327,SRR16648328,SRR16648329,SRR16648330,SRR16648331,SRR16648332,Unknown_BA854-01T0001,Unknown_BA854-01T0002,Unknown_BA854-01T0003 --rnaseq_sets_dirs=/scratch/project_2009273/RNAseq/ --prot_seq=./Viridiplantae.fa --eval_pseudo=pseudo.gff3 --AUGUSTUS_ab_initio --grass --gff3 --threads=$OMP_NUM_THREADS --busco_lineage=poales_odb10
#create a gff3 file
#cat braker.gtf | perl -ne 'if(m/\tAUGUSTUS\t/ or m/\tGeneMark\.hmm\t/ or m/\tGeneMark\.hmm3\t/ or m/\tgmst\t/) {print $_;}' | gtf2gff.pl --gff3 --out=braker.gff3
#script from EVM,
perl ./braker_GTF_to_EVM_GFF3.pl  braker.gtf > braker.evm_mod.gff3

#trinity, on mahti
cd /scratch/trinity
ln -s /scratch/RNAseq/merged_20.bam .
singularity exec -B /scratch/RNAseq:/scratch/RNAseq -B ${PWD}:${PWD} -e trinityrnaseq.v2.15.1.simg  Trinity --genome_guided_bam merged.bam --genome_guided_max_intron 10000 --max_memory 230G --CPU 128
#use array on puhti
cd /scratch/annotation/RNAseq
stringtie -o $name.out.gtf -G ../braker.gff3 -p 40 -l $name Aligned.out.bam.sort
awk '{ print "./"$1"/"$1".out.gtf" }' name.list > stringtiemerge.list
stringtie --merge -G braker.gff3 -o merge.gtf -l merged -p 40 stringtiemerge.list
perl gtf_genome_to_cdna_fasta.pl merge.gtf /scratch/gapclose/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta > transcripts.fasta
cat transcripts.fasta Trinity-GG.fasta > transcripts_all.fasta
#cyrus
conda activate ltr_retriever 
cd-hit-est -i transcripts_all.fasta -o transcripts_all.cdhit.fasta -n 10 -c 0.95 -d 0 -M 20000 -T 40
#mahti
cd /scratch/pasa
singularity exec -B ${PWD}:${PWD} pasapipeline.v2.5.3.simg seqclean transcripts_all.cdhit.fasta
ln -s /scratch/genome/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
module load sqlite
ORIGINALTRANSCRIPTOME=$(ls *.clean | sed 's/.clean//g')
CLEANTRANSCRIPTOME=$(ls *.clean)
CPU=$SLURM_CPUS_PER_TASK
REFERENCE_GENOME=/scratch/genome/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta
singularity exec -B ${PWD}:${PWD} -B /scratch/genome/:/scratch/genome/ pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c alignment_config.txt -C -R -g $REFERENCE_GENOME --ALIGNERS gmap,blat,minimap2 -t $CLEANTRANSCRIPTOME -T -u $ORIGINALTRANSCRIPTOME --CPU 32 &> alignment_assembly.log
singularity exec -B ${PWD}:${PWD} pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta PhragRNA.sqlite.assemblies.fasta  --pasa_transcripts_gff3 PhragRNA.sqlite.pasa_assemblies.gff3 &> transdecoder.log
#the braker final results contains too few genes, use the augustus and gnemarketp results before merging
perl /opt/EVidenceModeler/EvmUtils/misc/genemark_gtf2gff3.pl genemark.gtf > genemark_evm_mod.gff3 
#take the script from evidencemodlder/evmutils, augustus_GTF_to_EVM_GFF3.pl
/opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl augustus.hints.gtf > augustus_evm_mod.gff3

#mahti,install evidencemodler
cd /scratch/EVM
export SINGULARITY_TMPDIR=/scratch/EVM
export SINGULARITY_CACHEDIR=/scratch/EVM
unset XDG_RUNTIME_DIR
#wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/MISC/EVidenceModeler/EVidenceModeler.v2.1.0.simg
ln -s /scratch/genome/Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta .
cp /scratch/pasa/PhragRNA.sqlite.pasa_assemblies.gff3 .
cp /scratch/pasa/PhragRNA.sqlite.assemblies.fasta.transdecoder.genome.gff3 .
#cat genemark_evm_mod.gff3 braker.evm_mod.gff3 > gene_predictions.gff3 #NEW, do not use genemarkES for annotation any more, the IGV shows genemerkES is very inaccurate and cause fusing genes
#use the genmark and augustus gtf files produced by braker3
cat genemark_evm_mod.gff3 augustus_evm_mod.gff3> gene_predictions.gff3
cat Phalli_gemoma_evm_mod.gff3 Osativa_gemoma_evm_mod.gff3 Sviridis_gemoma_evm_mod.gff3 Sbicolor_gemoma_evm_mod.gff3  >  protein_alignments.gff3
cat PhragRNA.sqlite.pasa_assemblies.gff3 PhragRNA.sqlite.assemblies.fasta.transdecoder.genome.gff3 > transcriptome_alignments.gff3
##pay attention to the capital letters, use tab instead of space
echo "PROTEIN   GeMoMa   1" >> weights.txt
echo "ABINITIO_PREDICTION   GeneMark.hmm3   1" >> weights.txt
echo "ABINITIO_PREDICTION   AUGUSTUS   2" >> weights.txt
echo "ABINITIO_PREDICTION   gmst   1" >> weights.txt
echo "TRANSCRIPT   assembler-PhragRNA.sqlite   10" >> weights.txt
echo "OTHER_PREDICTION   transdecoder   7" >> weights.txt

#use cyrus
/opt/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcriptome_alignments.gff3 --repeats Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.out.gff --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
/opt/EVidenceModeler/EvmUtils/write_EVM_commands.pl --genome Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta --weights $PWD"/weights.txt" --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcriptome_alignments.gff3 --repeats Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.out.gff --output_file_name evm.out  --partitions partitions_list.out >  commands.list
linenumbers=$(wc -l commands.list | cut -d " " -f1) #Obtaining the number of lines in commands.list
splitnumber=$((linenumbers/40)) #I am dividing the line numbers by the number of CPUs given so that we can have almost equal number of lines in each file.
split -l $splitnumber commands.list #Carrying out the actual split
array1=($(ls x*)) #Initiating the loop for running EvM in parallel
for i in ${array1[@]}
do
        echo "Processing EVidenceModeler in parallel"
        /opt/EVidenceModeler/EvmUtils/execute_EVM_commands.pl $i | tee $i".log" &
done
wait
echo "Completed EVidenceModeler in parallel"
/opt/EVidenceModeler/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
/opt/EVidenceModeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out  --genome Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
/opt/EVidenceModeler/EvmUtils/gff3_file_to_proteins.pl EVM.all.gff3 Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta > EVM.all.faa

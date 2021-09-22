## MAKER annotaion for Erysiphe necator genome 

#### Build a transcriptome assembly

All SRR numbers for Erysiphe necator RNA data available on NCBI were written into a file called Ene_SRR_list.txt, one number per line.

File get_SRR_data.sh was created to include

```   
   #!/usr/bin/bash
   fasterq-dump $1
```

```ShellSession
cat Ene_SRR_list.txt | xargs -n 1 bash get_SRR_data.sh
```

Fastq-dump will pull the data, one by one for all accession numbers in your list, and turn each into a fastq at the same time.

```ShellSession
source /programs/HISAT2/hisat2.sh
hisat2-build -f Erysiphe_necator.fa Ene_index

hisat2 -p 40 -q -x Ene_index SRR13696274.fastq,SRR13696275.fastq,SRR13696276.fastq,SRR13696277.fastq,SRR13696278.fastq,SRR13696279.fastq,SRR13696280.fastq,SRR13696281.fastq,SRR13696282.fastq,SRR13696283.fastq,SRR13696284.fastq,SRR13696285.fastq,SRR13696286.fastq,SRR13696287.fastq,SRR13696288.fastq,SRR13696289.fastq,SRR13696290.fastq,SRR13696291.fastq,SRR13696292.fastq,SRR13696293.fastq,SRR13696294.fastq,SRR13696295.fastq,SRR13696296.fastq,SRR13696297.fastq,SRR13696298.fastq,SRR13696299.fastq,SRR13696300.fastq,SRR13696301.fastq,SRR13696302.fastq,SRR13696303.fastq,SRR13696304.fastq,SRR13696305.fastq,SRR13696306.fastq,SRR13696307.fastq,SRR13696308.fastq,SRR13696309.fastq,SRR13696310.fastq,SRR13696311.fastq,SRR13696312.fastq,SRR13696313.fastq,SRR13696314.fastq,SRR13696315.fastq,SRR13696316.fastq,SRR13696317.fastq,SRR13696318.fastq,SRR13696319.fastq,SRR13696320.fastq,SRR13696321.fastq --al Ene_alignedRNA -S Ene_alignment1
```

Replace spaces in sequence headers with underscore
```ShellSession
cat Ene_alignedRNA | perl -lane 's/\s/_/g; print;' > Ene_alignedRNA_V2.fastq
```

Run Trinity
```ShellSession 
export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib 
screen
$TRINITY_HOME/Trinity --seqType fq --single Ene_alignedRNA_V2.fastqs--SS_lib_type R --max_memory 160G --no_salmon --trimmomatic --CPU 40 --output ./trinity_out >& trinity.log &
```

The output file "Trinity_Ene.fasta" is to be supplied to "est=" in the MAKER control file.

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name One Erysiphe_necator.fa
RepeatModeler -pa 40 -database Ene -LTRStruct >& repeatmodeler.log
```
The output file "Ene-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

#### First `MAKER` annotation run

Set environment to run MAKER and create MAKER control files.

```ShellSession
cd /workdir
cp -r /programs/Augustus-3.3.2/config/ /workdir/
export PATH=/workdir/maker/bin:$PATH
export PATH=/programs/snap:$PATH
export AUGUSTUS_CONFIG_PATH=/workdir/config
export ZOE=/programs/snap/Zoe
export LD_LIBRARY_PATH=/programs/boost_1_62_0/lib
which maker          #The step is to confirm that you are using maker on /workdir
mkdir /workdir/tmp
cd nv232
maker -CTL
```

Edit maker_opts.ctl file

```
genome=/workdir/nv232/Erysiphe_necator.fa
organism_type=eukaryotic
est=/workdir/nv232/Trinity_Ene.fasta #output from trinity
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Ene-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model produced by @StefanKusch and moved to ./config/species
est2genome=1 #infer gene predictions directly from EST
min_contig=500
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base Ene_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../Ene_rnd1.maker.output/Ene_rnd1_master_datastore_index.log
maker2zff -l 50 -x 0.5 Ene_rnd1.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Ene . > ../snap_Ene_1.hmm
mv Ene_rnd1.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
maker_gff= Ene_rnd1.all.gff 
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
protein_pass=1 
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= #not running repeat masking again
repeat_protein= #not running repeat masking again
snaphmm=snap_Ene_1.hmm
augustus_species= #not running AUGUSTUS again
est2genome=0 # do not do EST evidence based gene model
min_contig=500
keep_preds=1
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Ene_rnd2 -fix_nucleotides -qq >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../Ene_rnd2.maker.output/Ene_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 Ene_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Ene . > ../snap_Ene_2.hmm
mv Ene_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

```
maker_gff=Ene_rnd2.all.gff
snaphmm=snap_Ene_2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Ene_rnd3 -fix_nucleotides -qq >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -d Ene_rnd3.maker.output/Ene_rnd3_master_datastore_index.log
fasta_merge -d Ene_rnd3.maker.output/Ene_rnd3_master_datastore_index.log
```

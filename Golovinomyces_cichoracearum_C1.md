## MAKER annotaion for Golovinomyces cichoracearum C1 genome 

#### Build a transcriptome assembly

```ShellSession
fasterq-dump --split-files SRR4048996.1
fasterq-dump --split-files SRR4048997.1
fasterq-dump --split-files SRR4048998.1
fasterq-dump --split-files SRR4048999.1

source /programs/HISAT2/hisat2.sh
hisat2-build -f GcichoracearumC1.genome.fa Gci_C1_index
hisat2 -p 40 -q -x Gci_C1_index -1 SRR4048996.1_1.fastq -2 SRR4048996.1_2.fastq --al-conc Gci_C1_alignedRNA_1 -S Gci_C1_alignment1
hisat2 -p 40 -q -x Gci_C1_index -1 SRR4048997.1_1.fastq -2 SRR4048997.1_2.fastq --al-conc Gci_C1_alignedRNA_2 -S Gci_C1_alignment2
hisat2 -p 40 -q -x Gci_C1_index -1 SRR4048998.1_1.fastq -2 SRR4048998.1_2.fastq --al-conc Gci_C1_alignedRNA_3 -S Gci_C1_alignment3
hisat2 -p 40 -q -x Gci_C1_index -1 SRR4048999.1_1.fastq -2 SRR4048999.1_2.fastq --al-conc Gci_C1_alignedRNA_4 -S Gci_C1_alignment4

export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib
/programs/trinityrnaseq-v2.10.0/Trinity --seqType fq --left Gci_C1_alignedRNA_1.1,Gci_C1_alignedRNA_2.1,Gci_C1_alignedRNA_3.1,Gci_C1_alignedRNA_4.1 /
--right Gci_C1_alignedRNA_1.2,Gci_C1_alignedRNA_2.2,Gci_C1_alignedRNA_3.2,Gci_C1_alignedRNA_4.2 --SS_lib_type RF --max_memory 256G --trimmomatic /
--CPU 40 --output ./trinity_out >& trinity.log &
```

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name Gci_C1 GcichoracearumC1.genome.fa
RepeatModeler -pa 40 -database Gci_C1 -LTRStruct >& repeatmodeler.log
```
The output file "Gci_C1-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

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
maker -CTL
```

Edit maker_opts.ctl file

```
genome=/workdir/nv232/GcichoracearumC1.genome.fa
organism_type=eukaryotic
est=/workdir/nv232/Trinity_Gci_C1.fasta #output from trinity
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Gci_C1-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model prduced by @StefanKusch and moved to ./config/species
est2genome=1 #infer gene predictions directly from EST
min_contig=500
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 64 maker -fix_nucleotides -base GciC1_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../GciC1_rnd1.genome.maker.output/GciC1_rnd1.genome_master_datastore_index.log
maker2zff -l 50 -x 0.5 GciC1_rnd1.genome.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl GciC1 . > ../snap_GciC1_1.hmm
mv GciC1_rnd1.genome.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
maker_gff= GciC1_rnd1.genome.all.gff 
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= #not running repeat masking again
repeat_protein= #not running repeat masking again
snaphmm=snap_GciC1_1.hmm
augustus_species= #not running AUGUSTUS again
est2genome=0 # do not do EST evidence based gene model
min_contig=500
keep_preds=1
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 64 maker -base GciC1_rnd2 -fix_nucleotides -qq >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../GciC1_rnd2.maker.output/GciC1_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 GciC1_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl GciC1 . > ../snap_GciC1_2.hmm
mv GciC1_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

```
maker_gff=GciC1_rnd2.all.gff
snaphmm=snap_GciC1_2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 64 maker -base GciC1_rnd3 -fix_nucleotides -qq >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -d GciC1_rnd3.maker.output/GciC1_rnd3_master_datastore_index.log
fasta_merge -d GciC1_rnd3.maker.output/GciC1_rnd3_master_datastore_index.log
```

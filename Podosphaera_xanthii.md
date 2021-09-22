## MAKER annotation for Podosphaera xanthii Wanju genome 

#### Build a transcriptome assembly

```ShellSession
fasterq-dump --split-files SRR8216439.1
fasterq-dump --split-files SRR8216440.1
fasterq-dump --split-files SRR8216441.1
fasterq-dump --split-files SRR8858150.1


source /programs/HISAT2/hisat2.sh
hisat2-build -f Podosphaera_xanthii_Wanju.fna PxaW_index
hisat2 -p 40 -q -x PxaW_index -1 SRR8216439.1_1.fastq -2 SRR8216439.1_2.fastq --al-conc PxaW_alignedRNA_1 -S PxaW_alignment1
hisat2 -p 40 -q -x PxaW_index -1 SRR8216440.1_1.fastq -2 SRR8216440.1_2.fastq --al-conc PxaW_alignedRNA_2 -S PxaW_alignment2
hisat2 -p 40 -q -x PxaW_index -1 SRR8216441.1_1.fastq -2 SRR8216441.1_2.fastq --al-conc PxaW_alignedRNA_3 -S PxaW_alignment3
hisat2 -p 40 -q -x PxaW_index -1 SRR8858150.1_1.fastq -2 SRR8858150.1_2.fastq --al-conc PxaW_alignedRNA_4 -S PxaW_alignment4
```

Replace spaces in sequence headers with underscore
```ShellSession
cat PxaW_alignedRNA_1.1 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_1.1_V2.fastq
cat PxaW_alignedRNA_1.2 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_1.2_V2.fastq
cat PxaW_alignedRNA_2.1 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_2.1_V2.fastq
cat PxaW_alignedRNA_2.2 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_2.2_V2.fastq
cat PxaW_alignedRNA_3.1 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_3.1_V2.fastq
cat PxaW_alignedRNA_3.2 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_3.2_V2.fastq
cat PxaW_alignedRNA_4.1 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_4.1_V2.fastq
cat PxaW_alignedRNA_4.2 | perl -lane 's/\s/_/g; print;' > PxaW_alignedRNA_4.2_V2.fastq
```

Run Trinity
```ShellSession
export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib 
screen
$TRINITY_HOME/Trinity --seqType fq --left PxaW_alignedRNA_1.1_V2.fastq,PxaW_alignedRNA_2.1_V2.fastq,PxaW_alignedRNA_3.1_V2.fastq,PxaW_alignedRNA_4.1_V2.fastq --right PxaW_alignedRNA_1.2_V2.fastq,PxaW_alignedRNA_2.2_V2.fastq,PxaW_alignedRNA_3.2_V2.fastq,PxaW_alignedRNA_4.2_V2.fastq --SS_lib_type RF --max_memory 160G --no_salmon --trimmomatic --CPU 64 --output ./trinity_out >& trinity.log &
```

The output file "Trinity.fasta" is to be supplied to "est=" in the MAKER control file.

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name pxw Podosphaera_xanthii_Wanju.fna
RepeatModeler -pa 40 -database pxw -LTRStruct >& repeatmodeler.log
```
The output file "pxw-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

#### First `MAKER` annotation run

Set environment to run MAKER and create MAKER control files.

```ShellSession
cd /workdir
cp -r /programs/Augustus-3.3.2/config/ /workdir/
cp -rH /programs/maker/ ./
cp -rH /programs/RepeatMasker  ./
export PATH=/workdir/$USER/maker/bin:/workdir/$USER/RepeatMasker:/programs/snap:$PATH
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
genome=/workdir/nv232/Podosphaera_xanthii_Wanju.fna
organism_type=eukaryotic
est=/workdir/nv232/Trinity_PxaW.fasta #output from trinity
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/pxw-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model produced by @StefanKusch and moved to ./config/species
est2genome=1 #infer gene predictions directly from EST
min_contig=500
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 64 maker -fix_nucleotides -base Pxw_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../Pxw_rnd1.maker.output/Pxw_rnd1_master_datastore_index.log
maker2zff -l 50 -x 0.5 Pxw_rnd1.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Pxw . > ../snap_Pxw_1.hmm
mv Pxw_rnd1.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
maker_gff= Pxw_rnd1.all.gff 
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
perotein_pass=1
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= #not running repeat masking again
repeat_protein= #not running repeat masking again
snaphmm=snap_Pxw_1.hmm
augustus_species= #not running AUGUSTUS again
est2genome=0 # do not do EST evidence based gene model
min_contig=500
keep_preds=1
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 64 maker -base Pxw_rnd2 -fix_nucleotides -qq >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../Pxw_rnd2.maker.output/Pxw_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 Pxw_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Pxw . > ../snap_Pxw_2.hmm
mv Pxw_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

```
maker_gff=Pxw_rnd2.all.gff
snaphmm=snap_Pxw_2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 64 maker -base Pxw_rnd3 -fix_nucleotides -qq >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -d Pxw_rnd3.maker.output/Pxw_rnd3_master_datastore_index.log
fasta_merge -d Pxw_rnd3.maker.output/Pxw_rnd3_master_datastore_index.log
```

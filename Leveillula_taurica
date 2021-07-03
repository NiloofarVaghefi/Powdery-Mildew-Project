## MAKER annotaion for Leveillula taurica genome 

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name Lta Leveillula_taurica_v1.fa
RepeatModeler -pa 40 -database Lta -LTRStruct >& repeatmodeler.log
```
The output file "Lta-families.fa" is to be supplied to "rmlib=" in the MAKER control file.


#### 'Trinity' transcriptome assembly

```ShellSession
/programs/bbmap-38.86/bbduk.sh in1=C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R1.fastq.gz in2=C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R2.fastq.gz out1=C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R1_BBduk.fastq.gz out2=C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R2_BBduk.fastq.gz ref=UniVec_Core ktrim=r k=21 mink=11 hdist=2 tpe tbo
/programs/bbmap-38.86/current/ jgi.BBDuk in1=C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R1.fastq.gz in2=C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R2.fastq.gz out1=C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R1_BBduk.fastq.gz out2=C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R2_BBduk.fastq.gz ref=UniVec_Core ktrim=r k=21 mink=11 hdist=2 tpe tbo 
export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib
screen
$TRINITY_HOME/Trinity --seqType fq --left C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R1_BBduk.fastq.gz,C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R1_BBduk.fastq.gz --right C1_HWVJKAFXY_TTGCCTAG-TAAGTGGT_R2_BBduk.fastq.gz,C2_HWVJKAFXY_CCATTCGA-CGGACAAC_R2_BBduk.fastq.gz --SS_lib_type RF --max_memory 160G --no_salmon --trimmomatic --CPU 40 --output ./Lta_trinity_out >& trinity.log &
```
 
 The output file "Trinity.fasta" is to be supplied to "rmlib=" in the MAKER control file.
 
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
genome=/workdir/nv232/Leveillula_taurica_v1.fa
est=/workdir/nv232/Trinity.fasta
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Lta-families.fa #organism specific repeat library output from RepeatModeler 
softmask = 1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model prduced by @StefanKusch
est2genome=1 #infer gene predictions directly from EST
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -base Lta_rnd1  >& log1 &
```

#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../Lta_rnd1.maker.output/Lta_rnd1_master_datastore_index.log
maker2zff -l 50 -x 0.5 Lta_rnd1.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Lta . > ../Lta1.hmm
mv Lta_rnd1.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
genome=/workdir/nv232/Leveillula_taurica_v1.fa
maker_gff=/workdir/nv232/Lta_rnd1.all.gff
est_pass=1 #use est alignment from round 1
protein_pass=1 #use protein alignment from round 1
rm_pass=1 # use repeats in the gff file
est= #remove
model_org= #remove
rmlib = #remove
repeat_protein= #Remove
snaphmm=/workdir/nv232/Lta1.hmm
est2genome=0
pred_stats=1 #report AED stats
keep_preds=1 # keep genes even without evidence support, set to 0 if no
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Lta_rnd2 >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../Lta_rnd2.maker.output/Lta_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 Lta_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Lta . > ../Lta2.hmm
mv Lta_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

``` 
maker_gff=Lta_rnd2.all.gff
snaphmm=Lta2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Lta_rnd3 >& log3 &
```

#### `MAKER` outputs
```ShellSession
gff3_merge -d Lta_rnd3.maker.output/Lta_rnd3_master_datastore_index.log>epi_rnd3.gff
fasta_merge -d Lta_rnd3.maker.output/Lta_rnd3_master_datastore_index.log
```

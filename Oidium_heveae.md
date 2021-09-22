## MAKER annotation for Oidium heveae genome 

#### Build a transcriptome assembly

```ShellSession
fasterq-dump --split-files SRR7716073.1
fasterq-dump --split-files SRR7716074.1
fasterq-dump --split-files SRR7716075.1
fasterq-dump --split-files SRR7716076.1

source /programs/HISAT2/hisat2.sh
hisat2-build -f Oidium_heveae.fa Ohe_index
hisat2 -p 64 -q -x Ohe_index -1 SRR7716073.1_1.fastq -2 SRR7716073.1_2.fastq --al-conc Ohe_alignedRNA_1 -S Ohe_alignment1
hisat2 -p 64 -q -x Ohe_index -1 SRR7716074.1_1.fastq -2 SRR7716074.1_2.fastq --al-conc Ohe_alignedRNA_2 -S Ohe_alignment2
hisat2 -p 64 -q -x Ohe_index -1 SRR7716075.1_1.fastq -2 SRR7716075.1_2.fastq --al-conc Ohe_alignedRNA_3 -S Ohe_alignment3
hisat2 -p 64 -q -x Ohe_index -1 SRR7716076.1_1.fastq -2 SRR7716076.1_2.fastq --al-conc Ohe_alignedRNA_4 -S Ohe_alignment4
```

Replace spaces in sequence headers with underscore
```ShellSession
cat Ohe_alignedRNA_1.1 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_1.1_V2.fastq
cat Ohe_alignedRNA_1.2 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_1.2_V2.fastq
cat Ohe_alignedRNA_2.1 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_2.1_V2.fastq
cat Ohe_alignedRNA_2.2 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_2.2_V2.fastq
cat Ohe_alignedRNA_3.1 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_3.1_V2.fastq
cat Ohe_alignedRNA_3.2 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_3.2_V2.fastq
cat Ohe_alignedRNA_4.1 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_4.1_V2.fastq
cat Ohe_alignedRNA_4.2 | perl -lane 's/\s/_/g; print;' > Ohe_alignedRNA_4.2_V2.fastq
```

Run Trinity
```ShellSession
export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib
screen
$TRINITY_HOME/Trinity --seqType fq --left Ohe_alignedRNA_1.1_V2.fastq,Ohe_alignedRNA_2.1_V2.fastq,Ohe_alignedRNA_3.1_V2.fastq,Ohe_alignedRNA_4.1_V2.fastq --right Ohe_alignedRNA_1.2_V2.fastq,Ohe_alignedRNA_2.2_V2.fastq,Ohe_alignedRNA_3.2_V2.fastq,Ohe_alignedRNA_4.2_V2.fastq --SS_lib_type RF --max_memory 160G --no_salmon --trimmomatic --CPU 40 --output ./trinity_out >& trinity.log &
```

The output file "Trinity.fasta" is to be supplied to "est=" in the MAKER control file.

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name Ohe Oidium_heveae.fa
RepeatModeler -pa 64 -database Ohe -LTRStruct >& repeatmodeler.log
```
The output file "Ohe-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

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
genome=/workdir/nv232/Oidium_heveae.fa
organism_type=eukaryotic
est=/workdir/nv232/Trinity_Ohe.fasta #output from trinity
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Ohe-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model prduced by @StefanKusch and moved to ./config/species
est2genome=1 #infer gene predictions directly from EST
min_contig=500
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base Ohe_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../Ohe_rnd1.maker.output/Ohe_rnd1_master_datastore_index.log
maker2zff -l 50 -x 0.5 Ohe_rnd1.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Ohe . > ../snap_Ohe1.hmm
mv Ohe_rnd1.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
maker_gff= Ohe_rnd1.all.gff 
est_pass=1 # use est alignment from round 1
protein_pass=1 # use protein alignment from round 1
rm_pass=1 # use repeats in the gff file
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= #not running repeat masking again
repeat_protein= #not running repeat masking again
snaphmm=snap_Ohe1.hmm
augustus_species= #not running AUGUSTUS again
est2genome=0 # do not do EST evidence based gene model
min_contig=500
keep_preds=1
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Ohe_rnd2 -fix_nucleotides -qq >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../Ohe_rnd2.maker.output/Ohe_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 Ohe_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Ohe . > ../snap_Ohe2.hmm
mv Ohe_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

```
maker_gff=Ohe_rnd2.all.gff
snaphmm=snap_Ohe2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Ohe_rnd3 -fix_nucleotides -qq >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -d Ohe_rnd3.maker.output/Ohe_rnd3_master_datastore_index.log
fasta_merge -d Ohe_rnd3.maker.output/Ohe_rnd3_master_datastore_index.log
```

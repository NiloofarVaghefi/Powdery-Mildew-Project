## MAKER annottaion for Erysiphe pisi genome 

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0:$PATH
BuildDatabase -name epi Erysiphe_pisi.genome.fa
RepeatModeler -pa 40 -database epi -LTRStruct >& repeatmodeler.log
```

#### First MAKER annotation run

```ShellSession
export PATH=/workdir/$USER/maker/bin:/workdir/$USER/RepeatMasker:/programs/snap:$PATH
export ZOE=/programs/snap/Zoe
export PATH=/programs/snap:$PATH
export LD_LIBRARY_PATH=/programs/boost_1_62_0/lib
cp -r /programs/Augustus-3.3.2/config/ /workdir/
export AUGUSTUS_CONFIG_PATH=/workdir/config
export PATH=/workdir/maker/bin:$PATH
which maker          #The step is to confirm that you are using maker on /workdir
mkdir /workdir/tmp
maker -CTL
```

Edit maker_opts.ctl file

```
genome=/workdir/nv232/Erysiphe_pisi.genome.fa
est=/workdir/nv232/Ep_GHEC01.1.fsa_nt
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Epi-families.fa #an organism specific repeat library in fasta format for RepeatMasker 
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering) 
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no 
TMP=/workdir/tmp	
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd1 -qq >& log &
```

#### Second MAKER annotation run - SNAP training

```
mkdir snap1
cd snap1
gff3_merge -d ../epi_rnd1.genome.maker.output/epi_rnd1.genome_master_datastore_index.log
maker2zff -l 50 -x 0.5 epi_rnd1.genome.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl epi . > ../snap_epi.hmm
mv epi_rnd1.genome.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1

Edit maker_opts.ctl file

```
maker_gff= epi_rnd1.genome.all.gff 
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
snaphmm=snap_epi.hmm
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= # not running repeat masking again
est2genome=0 # do not do EST evidence based gene model
pred_stats=1 #report AED stats
alt_splice=1 # 0: keep one isoform per gene; 1: identify splicing variants of the same gene
keep_preds=1 # keep genes even without evidence support, set to 0 if no
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd2 >& log2 &
```

#### Third MAKER annotation run - SNAP training

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../epi_rnd2.maker.output/epi_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 epi_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl epi . > ../snap_epi2.hmm
mv epi_rnd2.all.gff ..
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd2
```

Edit maker_opts.ctl file

```
maker_gff=epi_rnd2.all.gff
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
snaphmm=snap_spi2.hmm
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= # not running repeat masking again
est2genome=0 # do not do EST evidence based gene model
pred_stats=1 #report AED stats
alt_splice=1 # 0: keep one isoform per gene; 1: identify splicing variants of the same gene
keep_preds=1 # keep genes even without evidence support, set to 0 if no
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd3 >& log3 &
```

gff3_merge -n -d epi_rnd3.maker.output/epi_rnd3_master_datastore_index.log>epi_rnd3.noseq.gff
fasta_merge -d epi_rnd3.maker.output/epi_rnd3_master_datastore_index.log


rw-rw-r-- 1 nv232 nv232 322166388 Nov 20 22:55 log3
drwxrwxr-x 4 nv232 nv232       256 Nov 20 22:55 epi_rnd3.maker.output
-rw-rw-r-- 1 nv232 nv232         0 Nov 20 23:05 epi_rnd3.noseq.gff
-rw-rw-r-- 1 nv232 nv232 121960197 Nov 20 23:05 epi_rnd3.all.gff
-rw-rw-r-- 1 nv232 nv232  10692495 Nov 20 23:06 epi_rnd3.all.maker.snap_masked.transcripts.fasta
-rw-rw-r-- 1 nv232 nv232   4326207 Nov 20 23:07 epi_rnd3.all.maker.snap_masked.proteins.fasta
-rw-rw-r-- 1 nv232 nv232  12532791 Nov 20 23:07 epi_rnd3.all.maker.transcripts.fasta
-rw-rw-r-- 1 nv232 nv232   4362164 Nov 20 23:08 epi_rnd3.all.maker.proteins.fasta

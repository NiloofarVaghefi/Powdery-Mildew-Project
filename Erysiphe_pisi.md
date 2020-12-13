## MAKER annotaion for Erysiphe pisi genome 

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name epi Erysiphe_pisi.genome.fa
RepeatModeler -pa 40 -database epi -LTRStruct >& repeatmodeler.log
```
The output file "epi-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

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
genome=/workdir/nv232/Erysiphe_pisi.genome.fa
est=/workdir/nv232/Ep_GHEC01.1.fsa_nt #downloaded from NCBI TSA database
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Epi-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model prduced by @StefanKusch
est2genome=1 #infer gene predictions directly from EST
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

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
```

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
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd2 >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

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
snaphmm=snap_epi2.hmm
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= # not running repeat masking again
est2genome=0 # do not do EST evidence based gene model
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base epi_rnd3 >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -n -d epi_rnd3.maker.output/epi_rnd3_master_datastore_index.log>epi_rnd3.noseq.gff
fasta_merge -d epi_rnd3.maker.output/epi_rnd3_master_datastore_index.log
```

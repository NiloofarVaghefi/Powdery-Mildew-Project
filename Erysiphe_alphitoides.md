## MAKER annotaion for Erysiphe alphitoides genome 


#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name Eal Erysiphe_alphitoides.fasta
RepeatModeler -pa 40 -database Eal -LTRStruct >& repeatmodeler.log
```
The output file “Eal-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

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
export PATH=/workdir/$USER/maker/bin:/workdir/$USER/RepeatMasker:/programs/snap:$PATH
which maker          #The step is to confirm that you are using maker on /workdir
mkdir /workdir/nv232/tmp
cd nv232
maker -CTL
```

Edit maker_opts.ctl file

```
genome=/workdir/nv232/Erysiphe_alphitoides.fasta
organism_type=eukaryotic
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Eal-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model prduced by @StefanKusch and moved to ./config/species
min_contig=500
keep_preds=1
TMP=/workdir/nv232/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 64 maker -fix_nucleotides -base Eal_rnd1 -qq >& log &
```

#### `MAKER` outputs
```ShellSession
gff3_merge -d Eal_rnd1.maker.output/Eal_rnd1_master_datastore_index.log
fasta_merge -d Eal_rnd1.maker.output/Eal_rnd1_master_datastore_index.log
```

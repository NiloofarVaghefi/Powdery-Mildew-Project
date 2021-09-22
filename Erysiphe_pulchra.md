## MAKER annotation for Erysiphe pulchra genome 

#### Build a transcriptome assembly

All SRR numbers for Erysiphe pulchra RNA data available on NCBI were written into a file called SRR_Acc_List.txt, one number per line.

File get_SRR_data.sh was created to include

   #!/usr/bin/bash
   fasterq-dump $1

Fasterq-dump will pull the data, one by one for all accession numbers in your list, and turn each into a fastq at the same time. By default ,it will create paired end files if available.

```ShellSession
cat SRR_Acc_List5_Epu.txt | xargs -n 1 bash get_SRR_data.sh
cat *_1.fastq > left5.fastq
cat *_2.fastq > right5.fastq
rm SRR*
source /programs/HISAT2/hisat2.sh
hisat2-build -f Erysiphe_pulchra_TENN-F-071826.fna Epu_index
hisat2 -p 40 -q -x Epu_index -1 left5.fastq -2 right5.fastq --al-conc Epu_alignedRNA5 -S Epu_alignment5

cat SRR_Acc_List4_Epu.txt | xargs -n 1 bash get_SRR_data.sh
cat *_1.fastq > left4.fastq
cat *_2.fastq > right4.fastq
rm SRR*
source /programs/HISAT2/hisat2.sh
hisat2 -p 40 -q -x Epu_index -1 left4.fastq -2 right4.fastq --al-conc Epu_alignedRNA4 -S Epu_alignment4

cat SRR_Acc_List3_Epu.txt | xargs -n 1 bash get_SRR_data.sh
cat *_1.fastq > left3.fastq
cat *_2.fastq > right3.fastq
rm SRR*
source /programs/HISAT2/hisat2.sh
hisat2 -p 40 -q -x Epu_index -1 left3.fastq -2 right3.fastq --al-conc Epu_alignedRNA3 -S Epu_alignment3

cat SRR_Acc_List2_Epu.txt | xargs -n 1 bash get_SRR_data.sh
cat *_1.fastq > left2.fastq
cat *_2.fastq > right2.fastq
rm SRR*
source /programs/HISAT2/hisat2.sh
hisat2 -p 40 -q -x Epu_index -1 left2.fastq -2 right2.fastq --al-conc Epu_alignedRNA2 -S Epu_alignment2

cat SRR_Acc_List1_Epu.txt | xargs -n 1 bash get_SRR_data.sh
cat *_1.fastq > left1.fastq
cat *_2.fastq > right1.fastq
rm SRR*
source /programs/HISAT2/hisat2.sh
hisat2 -p 40 -q -x Epu_index -1 left1.fastq -2 right1.fastq --al-conc Epu_alignedRNA1 -S Epu_alignment1
```

Replace spaces in sequence headers with underscore
```ShellSession
cat Epu_alignedRNA5.2 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA5.2_V2.fastq
cat Epu_alignedRNA5.1 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA5.1_V2.fastq
cat Epu_alignedRNA4.2 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA4.2_V2.fastq
cat Epu_alignedRNA4.1 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA4.1_V2.fastq
cat Epu_alignedRNA3.1 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA3.1_V2.fastq
cat Epu_alignedRNA3.2 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA3.2_V2.fastq
cat Epu_alignedRNA2.1 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA2.1_V2.fastq
cat Epu_alignedRNA2.2 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA2.2_V2.fastq
cat Epu_alignedRNA1.1 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA1.1_V2.fastq
cat Epu_alignedRNA1.2 | perl -lane 's/\s/_/g; print;' > Epu_alignedRNA1.2_V2.fastq
```

Run Trinity
```ShellSession 
export PATH=/programs/jellyfish-2.2.7/bin:/programs/salmon-1.0.0/bin:$PATH
export TRINITY_HOME=/programs/trinityrnaseq-v2.10.0/
export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib 
screen
$TRINITY_HOME/Trinity --seqType fq --left Epu_alignedRNA5.1_V2.fastq,Epu_alignedRNA4.2_V2.fastq,Epu_alignedRNA3.1_V2.fastq,Epu_alignedRNA2.2_V2.fastq,Epu_alignedRNA1.1_V2.fastq --right Epu_alignedRNA5.2_V2.fastq,Epu_alignedRNA4.2_V2.fastq,Epu_alignedRNA3.2_V2.fastq,Epu_alignedRNA2.2_V2.fastq,Epu_alignedRNA1.2_V2.fastq --SS_lib_type RF --max_memory 160G --no_salmon --trimmomatic --CPU 40 --output ./trinity_out >& trinity.log &
```

The output file "Trinity_Epu.fasta" is to be supplied to "est=" in the MAKER control file.

#### Build a custom repeat database using `RepeatModeler`

```ShellSession
export PATH=/programs/RepeatModeler-2.0.1:$PATH
BuildDatabase -name Epu Erysiphe_pulchra_TENN-F-071826.fna
RepeatModeler -pa 40 -database Epu -LTRStruct >& repeatmodeler.log
```
The output file "Epu-families.fa" is to be supplied to "rmlib=" in the MAKER control file.

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
genome=/workdir/nv232/Erysiphe_pulchra_TENN-F-071826.fna
organism_type=eukaryotic
est=/workdir/nv232/Trinity_epu.fasta #output from trinity
model_org=simple #model organism for RepBase masking in RepeatMasker
rmlib=/workdir/nv232/Epu-families.fa #organism specific repeat library output from RepeatModeler 
softmask=1
augustus_species=Bgh_dh14_v4 #Augustus gene prediction species model produced by @StefanKusch and moved to ./config/species
est2genome=1 #infer gene predictions directly from EST
min_contig=500
TMP=/workdir/tmp
```

```ShellSession
screen
/usr/local/mpich/bin/mpiexec -n 40 maker -fix_nucleotides -base Epu_rnd1 -qq >& log &
```
 
#### `SNAP` training round 1 and second `MAKER` annotation run

```
mkdir snap1
cd snap1
gff3_merge -d ../Epu_rnd1.maker.output/Epu_rnd1_master_datastore_index.log
maker2zff -l 50 -x 0.5 Epu_rnd1.all.gff 
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Epu . > ../snap_Epu_1.hmm
mv Epu_rnd1.all.gff ../
cd ..
cp maker_opts.ctl  maker_opts.ctl_backup_rnd1
```

Edit maker_opts.ctl file

```
maker_gff= Epu_rnd1.all.gff 
est_pass=1 # use est alignment from round 1
rm_pass=1 # use repeats in the gff file
protein_pass=1 
est= # remove est file, do not run EST blast again
model_org= #remove repeat mask model, so not running RepeatModeler again
rmlib= #not running repeat masking again
repeat_protein= #not running repeat masking again
snaphmm=snap_Epu_1.hmm
augustus_species= #not running AUGUSTUS again
est2genome=0 # do not do EST evidence based gene model
min_contig=500
keep_preds=1
TMP=/workdir/tmp
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Epu_rnd2 -fix_nucleotides -qq >& log2 &
```

#### `SNAP` training round 2 and third `MAKER` annotation run 

```ShellSession
mkdir snap2
cd snap2
gff3_merge -d ../Epu_rnd2.maker.output/Epu_rnd2_master_datastore_index.log
maker2zff  -l 50 -x 0.5 Epu_rnd2.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Epu . > ../snap_Epu_2.hmm
mv Epu_rnd2.all.gff ..
cd ..
```

Edit maker_opts.ctl file

```
maker_gff=Epu_rnd2.all.gff
snaphmm=snap_Epu_2.hmm
```

```ShellSession
/usr/local/mpich/bin/mpiexec -n 40 maker -base Epu_rnd3 -fix_nucleotides -qq >& log3 &
```


#### `MAKER` outputs
```ShellSession
gff3_merge -d Epu_rnd3.maker.output/Epu_rnd3_master_datastore_index.log
fasta_merge -d Epu_rnd3.maker.output/Epu_rnd3_master_datastore_index.log
```

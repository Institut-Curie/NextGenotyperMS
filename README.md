# Table of Contents  
* [Description](#nextgenotyperms)  
* [Requirements](#requirements)
* [Usage](#usage)
* [Examples](#examples)
* [Output files](#output-files)
* [Limitations](#limitations)
* [Contact](#contact)
* [Citation](#citation)


# NextGenotyperMS
NextGenotyperMS is a tool written in python3 allowing to easily extract microsatellite genotypes from Next Generation Illumina reads. It has been developped in order to be able to easily compare sequencing results of microsatellites produced using a classical PCR approach against low temperature isothermal amplification using recombinase polymerase amplification (LT-RPA) ([Daunay et al. 2019](https://academic.oup.com/nar/article/47/21/e141/5570702])).

NextGenotyperMS can process fastq files or bams (ideally aligned with bwa) of paired samples and produce a high quality figure summarizing the distribution of the selected microsatellites. 

![Alt text](img/summary.png?raw=true "Microsatellite distribution")

# Requirements
NextGenotyperMS is distributed as a standalone [singularity](https://github.com/hpcng/singularity/releases) image so only singularity is required. It was tested on singularity version [3.6.1](https://github.com/hpcng/singularity/releases/tag/v3.6.1) but should be compatible with the higher versions of singularity as well.

# Usage
```
  -b BINDIR, --binDir=BINDIR
  -c CHRSIZEFILE, --chrSizeFile=CHRSIZEFILE
  -d DIRNAME, --dirName=DIRNAME
  -f FLANKSIZE, --flankSize=FLANKSIZE
  -F FIGFILENAME, --figFileName=FIGFILENAME
  --figSize=FIGSIZE     
  -m MINCOV, --minCov=MINCOV
  -M MINMAPQ, --minMapQ=MINMAPQ
  --maxKeyNb=MAXKEYNB   
  -n NBCPUS, --nbCpus=NBCPUS
  -p PLOTTARGETDIR, --plotTargetDir=PLOTTARGETDIR
  -q QUEUE, --queue=QUEUE
  -r REFFILE, --refFile=REFFILE
  -s SAMPLEFILE, --sampleFile=SAMPLEFILE
  -S SEQTOPLOTLIST, --seqToPlotList=SEQTOPLOTLIST
  --seqIdxToStartEndIdxDict=SEQIDXTOSTARTENDIDXDICT
  -t TARGETDIR, --targetDir=TARGETDIR
  -T TMPDIR, --tmpDir=TMPDIR
  --test=TEST           
  -u USEPERCENTINHIST, --usePercentInHist=USEPERCENTINHIST
  -x XLABELSIZE, --xLabelSize=XLABELSIZE
```


# Examples
```IMAGE=NextGenotyperMS.sif
TMP_DIR=/tmp
TARGET_DIR=/my/target/
NB_CPUS=4
REF_FILE=/usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/Sequence_MS_AmpSeq2.fa

singularity run $IMAGE nextGenotyper.py -d /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/ -r $REF_FILE -T $TMP_DIR -t $TARGET_DIR -n $NB_CPUS -s /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/Samples.txt -S MNRMS_HT17_145pb_T,MNRMS_NR24_128pb_T,MNRMS_CAT25_149pb_T,DNRMS_D2S123_227pb_CAxTA1CAy,QNRMS_REN_264pb_TCTG,QNRMS_HPRTII_304pb_TCTA --seqIdxToStartEndIdxDict "QNRMS_REN_264pb_TCTG:0;None,MNRMS_HT17_145pb_T:2;-3" -U 1 -F /data/tmp/vrenault/CEPH/testSet1_out/s.png --idvdToProcessList I2 --colorList green,yellow > /data/tmp/vrenault/CEPH/microsat/pipe.log
```

# Output files
* The distribution of each selected microsatellite in a single tab delimited file (the default name is ```$TARGET_DIR/\*_key.txt```) with two columns: the microsatellite length and the associated pattern as follows : 

```MNRMS_HT17_145pb_T
MS length       MS pattern
14      T14
15      T15
16      T16
17      T17
18      T18
19      T19
```

* one image in png format is created for each pair of individual and microsatellite showing the distribution of that microsatellite for the individual

* one figure combining the previous microsatellite distributions for all individuals is produced (thd default name is ```summary.tif```)


# Limitations
NextGenotyperMS will only work on microsatellites which size is < the sequencing read size.

# Contact
Victor RENAULT (victor.renault[AT]curie.fr)

# Citation

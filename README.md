# Table of Contents  
* [Description](#description)  
* [Requirements](#requirements)
* [Usage](#usage)
* [Examples](#examples)
* [Output files](#output-files)
* [Limitations](#limitations)
* [Contact](#contact)
* [Citation](#citation)


# Description
NextGenotyperMS is a tool written in python3 allowing to easily extract microsatellite genotypes from Next Generation Illumina reads. It has been developped in order to be able to easily compare sequencing results of microsatellites produced using a classical PCR approach against low temperature isothermal amplification using recombinase polymerase amplification (LT-RPA) ([Daunay et al. 2019](https://academic.oup.com/nar/article/47/21/e141/5570702])).

NextGenotyperMS can process fastq files or bams (ideally aligned with bwa) of paired (PCR/RPA) samples and produce a high quality figure summarizing the distribution of the selected microsatellites. 

![Alt text](img/summary.png?raw=true "Microsatellite distribution")

Briefly, each pair of fastq files goes through a round of QC using [fastqc (version 0.11.9)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and is then aligned using [bwa mem (version 0.7.17-r1188)](http://bio-bwa.sourceforge.net/). The generated bam or input bam is then *MarkDuplicated* (without duplicate removal) using [GATK (version 4.1.9.0)](https://github.com/broadinstitute/gatk/releases/tag/4.1.9.0) before genotypes are extracted by NextGenotyperMS using reads which position includes a given microsatellite listed in a given fasta reference sequence.

![Alt text](img/summary.png?raw=true "Microsatellite distribution")

# Requirements
NextGenotyperMS is distributed as a standalone [singularity](https://github.com/hpcng/singularity/releases) image so only singularity is required. It was tested on singularity version [3.6.1](https://github.com/hpcng/singularity/releases/tag/v3.6.1) but should be compatible with the higher versions of singularity as well.

NextGenotyperMS can run on a single machine and process all input data (fastq files or bams) or in a calculation cluster.
The minimum recommended RAM for the machine(s), independently of the analysis mode, is 8GB of RAM.

# Usage
```
singularity run  NextGenotyperMS.py -h
Usage: NextGenotyperMS.py [options]

Options:
  -h, --help            show this help message and exit
  -b BINDIR, --binDir=BINDIR
                        Folder containing all the necessary binaries. The
                        default value is "/usr/local/bin" and should not be
                        changed.
  -c CHRSIZEFILE, --chrSizeFile=CHRSIZEFILE
                        The file containing the length of the sequences in the
                        reference sequence (parameter "-r"). If left empty,
                        this file will be generated automatically.
  --colorList=COLORLIST
                        Set the matplotlib color for the PCR/RPA sample. The
                        default value is "red,royalblue" which sets PCR
                        samples to red and RPA samples to royalblue
  -d DIRNAME, --dirName=DIRNAME
                        The sample folder to process containing either a bam
                        file or paired fastq files. This option is only
                        relevant when bam/fastq files are processed
                        individually (when the option "-P local" is used)
                        typically on a calculation cluster
  -f FLANKSIZE, --flankSize=FLANKSIZE
                        the mininum number of bases that must be present at
                        each side of the microsatellite sequence extracted
                        from a read to be considered in the genotypes. The
                        default value is 3
  -F FIGFILENAME, --figFileName=FIGFILENAME
                        the output file for the microsatellite distribution
                        across all individuals. By default, the file will be
                        created in "TARGET_DIR/plots/summary.tif" where
                        TARGET_DIR is given by parameter "-t"
  --fastQcImageHeight=FASTQCIMAGEHEIGHT
                        Only relevant if your input data is fastq files. In
                        that case, this will set the height in pixels of
                        fastqc "per base sequence quality" images. The default
                        value is 200.
  --fastQcImageWidth=FASTQCIMAGEWIDTH
                        Only relevant if your input data is fastq files. In
                        that case, this will set the width in pixels of fastqc
                        "per base sequence quality" images. The default value
                        is 200.
  --figSize=FIGSIZE     comma separated value to set the width and height in
                        inches of the distribution figure. This is calculated
                        automatically by default
  -i IDVDTOPROCESSLIST, --idvdToProcessList=IDVDTOPROCESSLIST
                        comma-separated list of individuals to process. Set
                        only if you want to process only specific individuals
                        listed in the file given using the parameter "-s"
  -I IDVDTOEXCLUDELIST, --idvdToExcludeList=IDVDTOEXCLUDELIST
                        comma-separated list of individuals to exclude. Set
                        only if you want to exclude specific individuals
                        listed in the file given using the parameter "-s"
  -M MINMAPQ, --minMapQ=MINMAPQ
                        The minimum mapq score for a read to be considered in
                        the microsatellite genotypes. The default value is 20.
  -n NBCPUS, --nbCpus=NBCPUS
                        the maximum number of cores to use for one sample. If
                        you run NextGenotyperMS on a single machine, the value
                        should be set to the square root of the total nb of
                        cores of that machine rounded down to the nearest
                        integer. The default value is 1.
  -N FASTQCNBIMGPERLINE, --fastQcNbImgPerLine=FASTQCNBIMGPERLINE
                        Only relevant if your input data is fastq files. In
                        that case, this will set the number of sample fastqc
                        "per base sequence quality" images per line. The
                        default value is 5.
  -p PLOTTARGETDIR, --plotTargetDir=PLOTTARGETDIR
                        set the output folder for individual microsatellite
                        distributions and various summary files. By default,
                        this folder is set to TARGET_DIR/plots where
                        TARGET_DIR is given by parameter "-t"
  -P PROGNAME, --progName=PROGNAME
                        Leave empty if you plan to process all your data from
                        one single machine. Set to "prepareIdxForRef" in order
                        to create all the related index files for the
                        reference sequence given in parameter "-r". Set to
                        "local" when processing one sample at a time (useful)
                        for parallelization on a calculation cluster
  -r REFFILE, --refFile=REFFILE
                        The reference sequence in fasta format listing the
                        microsatellites of interest. The name of each
                        microsatellite sequence should finish with the
                        expected sequence pattern preceded by the character
                        "_". The supported patterns are a single repeated
                        sequence like "A", "TCTG" for example (eg
                        MNRMS_BAT26_195pb_A, QNRMS_REN_264pb_TCTG), two
                        patterns P and Q indicated as "PxQy" (eg
                        TNRMS_D12ATA63_247pb_TTGxTTAy. Please note that "x"
                        and "y" should be in lower case while the sequence
                        pattern should be in upper case) and microsatellites
                        with sequence patterns "PxQ1Ry" (eg
                        DNRMS_D2S123_227pb_CAxTA1CAy,
                        DNRMS_D18S61_157pb_CAxCG1CAy). An example can be found
                        in /usr/local/code3/curie/testNextGenotyperMS/Sequence
                        _MS_AmpSeq2.fa.
  -s SAMPLEFILE, --sampleFile=SAMPLEFILE
                        A tab-delimited file with at least the following 3
                        columns "SampleName" (should not contain any "_" nor
                        "."), "IdvdName" and "ExpType" (should be either "PCR"
                        or "RPA"). If the input files are in fastq format,
                        another column "FastqName" is necessary (fastq files
                        are expected to be named "{F}_R{N}_*.fastq.gz" or
                        "{F}.R{N}.*fastq.gz" where N is 1 or 2 and F should be
                        the value indicated in the column "FastqName"). An
                        example can be found in /usr/local/code3/curie/testNex
                        tGenotyperMS/bams/Samples.txt
  -S SEQTOPLOTLIST, --seqToPlotList=SEQTOPLOTLIST
                        comma-separated list of sequences listed in the
                        reference sequence (parameter "-r") to plot. The
                        default behaviour is to show all the sequences listed
                        in the reference sequence.
  --seqFileAlias=SEQFILEALIAS
                        tab-delimited file with two columns "seqName" (the
                        sequence name extracted from the reference sequence
                        file) and "alias" the sequence name to display
  -t TARGETDIR, --targetDir=TARGETDIR
                        the target folder where the analysis files should be
                        written to
  -T TMPDIR, --tmpDir=TMPDIR
                        the folder where temporary files should be written to
  -U USECUSTOMSEQNAME, --useCustomSeqName=USECUSTOMSEQNAME
                        if set to 1 the displayed sequence name will be the
                        second element extracted from the sequence name
                        splitted by "_". By default, the full sequence name is
                        shown
```


# Examples
## Single machine mode
As mentioned in the [usage](#usage) section, the parameter *-n* should be set to the square root of the total number of cores of your machine rounded down to the nearest integer. **Caution** : if you are not sure, leave it to 1 and **do not overestimate** that number as your machine will freeze if that number is too high
```IMAGE=NextGenotyperMS.sif
TMP_DIR=/tmp
TARGET_DIR=/my/target/
NB_CPUS=4 # in this example, my machine has 16 cores 
REF_FILE=/usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/Sequence_MS_AmpSeq2.fa
LOG_FILE=/my/log.txt

singularity run $IMAGE nextGenotyper.py -d /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/ -r $REF_FILE -T $TMP_DIR -t $TARGET_DIR -n $NB_CPUS -s /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/Samples.txt -U 1 > $LOG_FILE

singularity run $IMAGE nextGenotyper.py -d /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/ -r $REF_FILE -T $TMP_DIR -t $TARGET_DIR -n $NB_CPUS -s /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/Samples.txt -S MNRMS_HT17_145pb_T,MNRMS_NR24_128pb_T,MNRMS_CAT25_149pb_T,DNRMS_D2S123_227pb_CAxTA1CAy,QNRMS_REN_264pb_TCTG,QNRMS_HPRTII_304pb_TCTA -U 1 -F $TARGET_DIR/customSummary.png --idvdToProcessList I2 --colorList green,yellow > $LOG_FILE
```

## Calculation cluster
### parallelizable steps
```IMAGE=NextGenotyperMS.sif
TMP_DIR=/tmp
TARGET_DIR=/my/target/
NB_CPUS=4 # the number of available cores on the machine
REF_FILE=/usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/Sequence_MS_AmpSeq2.fa
LOG_FILE=/my/log.txt

singularity run $IMAGE nextGenotyper.py -d /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/ -r $REF_FILE -T $TMP_DIR -t $TARGET_DIR -n $NB_CPUS -s /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/Samples.txt -U 1 > $LOG_FILE
```


### summary step
Once all the parallelizable steps are finished, a final step is necessary to combine all the individual results.
```IMAGE=NextGenotyperMS.sif
TMP_DIR=/tmp
TARGET_DIR=/my/target/
NB_CPUS=4 # the number of available cores on the machine
REF_FILE=/usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/Sequence_MS_AmpSeq2.fa
LOG_FILE=/my/log.txt

singularity run $IMAGE nextGenotyper.py -d /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/ -r $REF_FILE -T $TMP_DIR -t $TARGET_DIR -s /usr/local/NextGenotyperMS/code3/curie/testNextGenotyperMS/testSet1/Samples.txt -U 1 > $LOG_FILE
```


# Output files
All result files are stored in the folder given to the parameter *-p*. If not specified, the default output folder is *TARGET_DIR/plots* where *TARGET_DIR* is given by parameter *-t*.

* one image in png format, which will be named *SEQNAME_INDIVIDUAL.png*, is created for each pair of individual and microsatellite showing the distribution of that microsatellite for the individual. A tab-delimited file with a *.txt* extension is also created with the same information.

* one figure combining the previous microsatellite distributions for all individuals is produced (the default name is ```summary.tif```). A tab-delimited file named by default ```summary_key.txt``` is also created with the correspondence between the microsatellite length shown in the figure and the related microsatellite sequence pattern.

* a tab-delimited file named ```mappinqQc.txt``` containing basic alignment metrics is available.

* if NextGenotyperMS is run from fastq files, a summary file named ```fastQcSummary.html``` is also created showing all the "per base sequence quality" images produced by fastqc.



# Limitations
* NextGenotyperMS will only work on microsatellites which size is < the sequencing read size.
* As described in the [usage](#usage) section, only microsatellites with the following patterns are handled:
** single repeated patterns like "A", "TCTG"
** two repeated patterns *P* and *Q* indicated as *PxQy* ie *P* and *Q* are respectively repeated *x* and *y* times
** microsatellites with sequence patterns *PxQ1Ry* ie *P*, *Q* and *R* are respectively repeated *x*, 1 and *y* times

# Contact
Victor RENAULT (victor.renault[AT]curie.fr)

# Citation

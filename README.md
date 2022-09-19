# xenocell_nextflow

Nextflow wrapper to execute XenoCell, a tool to separate graft from host reads in single-cell RNA-seq xenograft experiments (https://gitlab.com/XenoCell/XenoCell, https://doi.org/10.1186/s12920-021-00872). Not affiliated with the original authors of XenoCell.

Currently supports up to two lanes per fastq file.

**Command line options:**

```
--index: Combined Xenome reference genome index [currently needs to be pre-generated before running the workflow]
--outdir: Desired output directory
```

Method parameters (see documentation of the original tool):
```
--barcode_start: [Default: 1]
--barcode_length: = [Default: 16]
--lower_threshold_host: = [Default: 0.9]
--upper_threshold_graft: = [Default: 0.1]
--threads = [Default: 16]
--memory = [Default: 60]
--compression_level = [Default: 1]
```

Any other method parameters may be changed from defaults by modifying the code of xenocell.nf.
```
--fastqdir: The base directory storing sample fastq files. 
```

Currently expects the latter to be contained within subdirectories named after each sample within this base directory. Reads are expected to be named in a way that ends with a pattern like "_L001_R1_001.fastq.gz". To change the latter, modify the parameters below with any alternative pattern:
```
--fastq_lane_1_read_1: [Default: "${params.fastqdir}/*/*_L001_R1_001.fastq.gz"]
--fastq_lane_1_read_2: [Default: "${params.fastqdir}/*/*_L001_R2_001.fastq.gz"]
--fastq_lane_2_read_1: [Default: "${params.fastqdir}/*/*_L002_R1_001.fastq.gz"]
--fastq_lane_2_read_2: [Default: "${params.fastqdir}/*/*_L002_R2_001.fastq.gz"]
```

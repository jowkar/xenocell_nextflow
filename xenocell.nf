#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.threads = 1 // Might be safest, since Xenome has a known threading issue
params.memory = 24 // https://github.com/data61/gossamer/issues/14
params.compression_level = 1

params.graft_fasta = "/data/local/reference/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
params.host_fasta = "/data/local/reference/cellranger/refdata-gex-mm10-2020-A/fasta/genome.fa"
params.index = "/data/local/reference/cellranger/xenome/graft_GRCh38_host_mm10"
params.outdir = "./results"

params.barcode_start = 1
params.barcode_length = 16
params.lower_threshold_host = 0.9
params.upper_threshold_graft = 0.1

params.fastqdir = "/data/proj/um_ss/Investigations/data/raw/GWA-JN-605_PDX_expr/outs/fastq_path/GWA-JN-605"

params.fastq_lane_1_read_1 = "${params.fastqdir}/*/*_L001_R1_001.fastq.gz"
params.fastq_lane_1_read_2 = "${params.fastqdir}/*/*_L001_R2_001.fastq.gz"
params.fastq_lane_2_read_1 = "${params.fastqdir}/*/*_L002_R1_001.fastq.gz"
params.fastq_lane_2_read_2 = "${params.fastqdir}/*/*_L002_R2_001.fastq.gz"

process GENERATE_INDEX {

    input:
        past(graft_fasta)
        path(host_fasta)

    output:
        path(index), emit: index_ch

    script:
    """    
    xenocell.py generate_index \
        --host "${host_fasta}" \
        --graft "${graft_fasta}" \
        --output "${index}" \
        --threads 1 \
        --memory 24
    """
}

process MERGE_READS {

    input:
        tuple val(sample), val(lane_1_read_1), val(lane_2_read_1), val(lane_1_read_2), val(lane_2_read_2)
    
    output:
        tuple val(sample), path("*_R1.fastq.gz"), path("*_R2.fastq.gz"), emit: merge_reads_ch

    script:
    """
    cat ${lane_1_read_1} ${lane_2_read_1} > ${sample}_R1.fastq.gz
    cat ${lane_1_read_2} ${lane_2_read_2} > ${sample}_R2.fastq.gz
    """
}

process CLASSIFY_READS {

    // maxForks 1

    publishDir "${params.outdir}/xenocell/", mode: 'symlink'

    input:
        tuple val(sample), val(fastq_1), val(fastq_2)
        path(index)

    output:
        path("${sample}"), emit: classify_reads_ch

    script:
    """
    output="${sample}"
    if [ ! -d \$output ];
    then
        mkdir \$output
    fi

    xenocell.py classify_reads \
        --transcript "${fastq_2}" \
        --barcode "${fastq_1}" \
        --barcode_start "${params.barcode_start}" \
        --barcode_length "${params.barcode_length}" \
        --index "${params.index}" \
        --output "./" \
        --threads "${params.threads}" \
        --memory "${params.memory}" \
        --compression_level "${params.compression_level}"

    mv *.* \$output/
    """
}

process EXTRACT_GRAFT_BARCODES {
    
    // maxForks 1

    input:
        path(sample_output_dir)

    output:
        path("${sample_output_dir}/graft/*")

    script:
    """
    xenocell.py extract_cellular_barcodes \
        --input "${sample_output_dir}" \
        --barcode_start "${params.barcode_start}" \
        --barcode_length "${params.barcode_length}" \
        --subset_name graft \
        --lower_threshold 0.0 \
        --upper_threshold "${params.upper_threshold_graft}" \
        --threads "${params.threads}" \
        --compression_level "${params.compression_level}"
    """
}

process EXTRACT_HOST_BARCODES {
    
    // maxForks 1

    input:
        path(sample_output_dir)

    output:
        path("${sample_output_dir}/host/*")
    
    script:
    """
    xenocell.py extract_cellular_barcodes \
        --input "${sample_output_dir}" \
        --barcode_start "${params.barcode_start}" \
        --barcode_length "${params.barcode_length}" \
        --subset_name host \
        --lower_threshold "${params.lower_threshold_host}" \
        --upper_threshold 1.0 \
        --threads "${params.threads}" \
        --compression_level "${params.compression_level}"
    """
}

workflow {

    if (!params.index){
        GENERATE_INDEX(params.graft_fasta,params.host_fasta)
        index = GENERATE_INDEX.out.index_ch
    } else {
        index = Channel.fromPath(params.index)
    }

    //input_reads_ch = Channel.fromFilePairs(params.fastq, flat:true, size:4).view()
    input_reads_ch_lane_1_read_1 = Channel.fromFilePairs(params.fastq_lane_1_read_1, flat:true, size:1)
    input_reads_ch_lane_1_read_2 = Channel.fromFilePairs(params.fastq_lane_1_read_2, flat:true, size:1)
    input_reads_ch_lane_2_read_1 = Channel.fromFilePairs(params.fastq_lane_2_read_1, flat:true, size:1)
    input_reads_ch_lane_2_read_2 = Channel.fromFilePairs(params.fastq_lane_2_read_2, flat:true, size:1)

    input_reads_ch_read_1_joined = input_reads_ch_lane_1_read_1.join(input_reads_ch_lane_2_read_1)
    input_reads_ch_read_2_joined = input_reads_ch_lane_1_read_2.join(input_reads_ch_lane_2_read_2)

    input_reads_ch = input_reads_ch_read_1_joined.join(input_reads_ch_read_2_joined)

    MERGE_READS(
        input_reads_ch
    )

    CLASSIFY_READS(
         MERGE_READS.out.merge_reads_ch,
         index
    )

    EXTRACT_GRAFT_BARCODES(
        CLASSIFY_READS.out.classify_reads_ch
    )

    EXTRACT_HOST_BARCODES(
        CLASSIFY_READS.out.classify_reads_ch
    )
}
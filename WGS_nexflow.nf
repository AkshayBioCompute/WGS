// main.nf
params.ref = "/home/akshay/Akshay/WGS/hg38.fa.gz"
params.known_sites_vcf = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf"
params.known_sites_idx = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
params.aligned_reads = "/home/akshay/Akshay/WGS/aligned_reads"
params.reads = "/home/akshay/Akshay/WGS/reads"
params.results = "/home/akshay/Akshay/WGS/results"
params.data = "/home/akshay/Akshay/WGS/data"

process fastqc {
    input:
    path fastq from file(params.reads + "/SRR741411_*.filt.fastq.gz")

    output:
    path "${params.reads}/fastqc"

    script:
    """
    fastqc $fastq -o ${params.reads}
    """
}

process bwa_index {
    input:
    path ref

    output:
    path "${ref}.bwt"

    script:
    """
    bwa index $ref
    """
}

process bwa_mem {
    input:
    path ref, path r1, path r2

    output:
    path "${params.aligned_reads}/SRR741411.paired.sam"

    script:
    """
    bwa mem -t 4 -R '@RG\\tID:SRR741411\\tPL:ILLUMINA\\tSM:SRR741411' $ref $r1 $r2 > ${params.aligned_reads}/SRR741411.paired.sam
    """
}

process mark_duplicates {
    input:
    path samfile

    output:
    path "${params.aligned_reads}/SRR741411_sorted_dedup_reads.bam"

    script:
    """
    gatk MarkDuplicatesSpark -I $samfile -O ${params.aligned_reads}/SRR741411_sorted_dedup_reads.bam
    """
}

process base_recalibration {
    input:
    path bam, path ref, path known_sites

    output:
    path "${params.data}/recal_data.table"

    script:
    """
    gatk BaseRecalibrator -I $bam -R $ref --known-sites $known_sites -O ${params.data}/recal_data.table
    """
}

process apply_bqsr {
    input:
    path bam, path ref, path recal_table

    output:
    path "${params.aligned_reads}/SRR741411_sorted_dedup_bqsr_reads.bam"

    script:
    """
    gatk ApplyBQSR -I $bam -R $ref --bqsr-recal-file $recal_table -O ${params.aligned_reads}/SRR741411_sorted_dedup_bqsr_reads.bam
    """
}


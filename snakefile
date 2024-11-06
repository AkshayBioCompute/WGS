# Snakefile
ref = "/home/akshay/Akshay/WGS/hg38.fa.gz"
known_sites_vcf = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites_idx = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
aligned_reads = "/home/akshay/Akshay/WGS/aligned_reads"
reads = "/home/akshay/Akshay/WGS/reads"
results = "/home/akshay/Akshay/WGS/results"
data = "/home/akshay/Akshay/WGS/data"

rule all:
    input: 
        f"{results}/raw_snps.vcf", 
        f"{results}/raw_indels.vcf"

rule fastqc:
    input: 
        fastq="{reads}/{sample}_1.filt.fastq.gz"
    output: 
        qc_report="{reads}/{sample}_1_fastqc.zip"
    shell:
        "fastqc {input.fastq} -o {reads}"

rule bwa_index:
    input:
        ref
    output:
        ref + ".bwt"
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        ref=ref,
        r1=f"{reads}/SRR741411_1.filt.fastq.gz",
        r2=f"{reads}/SRR741411_2.filt.fastq.gz"
    output:
        f"{aligned_reads}/SRR741411.paired.sam"
    shell:
        "bwa mem -t 4 -R '@RG\\tID:SRR741411\\tPL:ILLUMINA\\tSM:SRR741411' {input.ref} {input.r1} {input.r2} > {output}"

rule mark_duplicates:
    input:
        f"{aligned_reads}/SRR741411.paired.sam"
    output:
        f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam"
    shell:
        "gatk MarkDuplicatesSpark -I {input} -O {output}"

rule base_recalibration:
    input:
        bam=f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam",
        ref=ref,
        known_sites=known_sites_vcf
    output:
        table=f"{data}/recal_data.table"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output.table}"
        
rule apply_bqsr:
    input:
        bam=f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam",
        ref=ref,
        table=f"{data}/recal_data.table"
    output:
        f"{aligned_reads}/SRR741411_sorted_dedup_bqsr_reads.bam"
    shell:
        "gatk ApplyBQSR -I {input.bam} -R {input.ref} --bqsr-recal-file {input.table} -O {output}"


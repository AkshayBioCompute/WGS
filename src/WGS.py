# src/WGS.py

import subprocess

# Define paths
ref = "/home/akshay/Akshay/WGS/hg38.fa.gz"
known_sites_vcf = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites_idx = "/home/akshay/Akshay/WGS/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
aligned_reads = "/home/akshay/Akshay/WGS/aligned_reads"
reads = "/home/akshay/Akshay/WGS/reads"
results = "/home/akshay/Akshay/WGS/results"
data = "/home/akshay/Akshay/WGS/data"

def run_fastqc():
    subprocess.run(["fastqc", f"{reads}/SRR741411_1.filt.fastq.gz", "-o", reads])
    subprocess.run(["fastqc", f"{reads}/SRR741411_2.filt.fastq.gz", "-o", reads])

def run_bwa_mem():
    subprocess.run(["bwa", "index", ref])
    subprocess.run([
        "bwa", "mem", "-t", "4", "-R", "@RG\\tID:SRR741411\\tPL:ILLUMINA\\tSM:SRR741411",
        ref, f"{reads}/SRR741411_1.filt.fastq.gz", f"{reads}/SRR741411_2.filt.fastq.gz",
        "-o", f"{aligned_reads}/SRR741411.paired.sam"
    ])

def mark_duplicates():
    subprocess.run([
        "gatk", "MarkDuplicatesSpark",
        "-I", f"{aligned_reads}/SRR741411.paired.sam",
        "-O", f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam"
    ])

def base_quality_recalibration():
    subprocess.run([
        "gatk", "BaseRecalibrator",
        "-I", f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam",
        "-R", ref, "--known-sites", known_sites_vcf,
        "-O", f"{data}/recal_data.table"
    ])
    subprocess.run([
        "gatk", "ApplyBQSR",
        "-I", f"{aligned_reads}/SRR741411_sorted_dedup_reads.bam",
        "-R", ref, "--bqsr-recal-file", f"{data}/recal_data.table",
        "-O", f"{aligned_reads}/SRR741411_sorted_dedup_bqsr_reads.bam"
    ])

def main():
    run_fastqc()
    run_bwa_mem()
    mark_duplicates()
    base_quality_recalibration()

if __name__ == "__main__":
    main()


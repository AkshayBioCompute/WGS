// WGS_Processing.wdl
workflow WGS_Processing {
    input {
        File ref                       # Path to reference genome, hg38.fa.gz
        File known_sites_vcf           # Path to known sites VCF, Homo_sapiens_assembly38.dbsnp138.vcf
        File known_sites_idx           # Index file for the known sites VCF, Homo_sapiens_assembly38.dbsnp138.vcf.idx
        File r1_fastq                  # Path to read 1 FASTQ, SRR741411_1.filt.fastq.gz
        File r2_fastq                  # Path to read 2 FASTQ, SRR741411_2.filt.fastq.gz
    }

    # Step 1: FastQC
    call FastQC as fastqc_r1 {
        input:
            fastq = r1_fastq
    }
    
    call FastQC as fastqc_r2 {
        input:
            fastq = r2_fastq
    }

    # Step 2: BWA Index
    call BWA_Index {
        input:
            ref = ref
    }

    # Step 3: BWA MEM
    call BWA_MEM {
        input:
            ref = ref,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq
    }

    # Step 4: Mark Duplicates
    call MarkDuplicates {
        input:
            sam = BWA_MEM.output_sam
    }

    # Step 5: Base Recalibration
    call BaseRecalibrator {
        input:
            bam = MarkDuplicates.output_bam,
            ref = ref,
            known_sites_vcf = known_sites_vcf
    }

    call ApplyBQSR {
        input:
            bam = MarkDuplicates.output_bam,
            ref = ref,
            recal_table = BaseRecalibrator.recal_table
    }

    # Step 6: Collect Alignment Metrics
    call CollectAlignmentMetrics {
        input:
            bam = ApplyBQSR.output_bam,
            ref = ref
    }

    # Step 7: Haplotype Caller
    call HaplotypeCaller {
        input:
            bam = ApplyBQSR.output_bam,
            ref = ref
    }

    # Step 8: Select Variants (SNPs)
    call SelectVariantsSNPs {
        input:
            ref = ref,
            vcf = HaplotypeCaller.raw_variants_vcf
    }

    # Step 8: Select Variants (Indels)
    call SelectVariantsIndels {
        input:
            ref = ref,
            vcf = HaplotypeCaller.raw_variants_vcf
    }
}

# Task definitions for each step of the workflow
task FastQC {
    input {
        File fastq
    }
    command {
        fastqc ${fastq} -o .
    }
    output {
        File qc_report = glob("*.zip")[0]
    }
    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
    }
}

task BWA_Index {
    input {
        File ref
    }
    command {
        bwa index ${ref}
    }
    output {
        Array[File] index_files = glob("${ref}*")
    }
    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
    }
}

task BWA_MEM {
    input {
        File ref
        File r1_fastq
        File r2_fastq
    }
    command {
        bwa mem -t 4 -R "@RG\\tID:SRR741411\\tPL:ILLUMINA\\tSM:SRR741411" ${ref} ${r1_fastq} ${r2_fastq} > SRR741411.paired.sam
    }
    output {
        File output_sam = "SRR741411.paired.sam"
    }
    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
    }
}

task MarkDuplicates {
    input {
        File sam
    }
    command {
        gatk MarkDuplicatesSpark -I ${sam} -O SRR741411_sorted_dedup_reads.bam
    }
    output {
        File output_bam = "SRR741411_sorted_dedup_reads.bam"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task BaseRecalibrator {
    input {
        File bam
        File ref
        File known_sites_vcf
    }
    command {
        gatk BaseRecalibrator -I ${bam} -R ${ref} --known-sites ${known_sites_vcf} -O recal_data.table
    }
    output {
        File recal_table = "recal_data.table"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task ApplyBQSR {
    input {
        File bam
        File ref
        File recal_table
    }
    command {
        gatk ApplyBQSR -I ${bam} -R ${ref} --bqsr-recal-file ${recal_table} -O SRR741411_sorted_dedup_bqsr_reads.bam
    }
    output {
        File output_bam = "SRR741411_sorted_dedup_bqsr_reads.bam"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task CollectAlignmentMetrics {
    input {
        File bam
        File ref
    }
    command {
        gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${bam} -O alignment_metrics.txt
    }
    output {
        File alignment_metrics = "alignment_metrics.txt"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task HaplotypeCaller {
    input {
        File bam
        File ref
    }
    command {
        gatk HaplotypeCaller -R ${ref} -I ${bam} -O raw_variants.vcf
    }
    output {
        File raw_variants_vcf = "raw_variants.vcf"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task SelectVariantsSNPs {
    input {
        File ref
        File vcf
    }
    command {
        gatk SelectVariants -R ${ref} -V ${vcf} --select-type SNP -O raw_snps.vcf
    }
    output {
        File raw_snps_vcf = "raw_snps.vcf"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}

task SelectVariantsIndels {
    input {
        File ref
        File vcf
    }
    command {
        gatk SelectVariants -R ${ref} -V ${vcf} --select-type INDEL -O raw_indels.vcf
    }
    output {
        File raw_indels_vcf = "raw_indels.vcf"
    }
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
    }
}


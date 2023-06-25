rule fastqc_raw_data:
    input:
        expand("/mnt/groupMansuy/kerem/tasks/longrna/rawdatas/HSC_{sample}.fastq.gz", sample=glob_wildcards("/mnt/groupMansuy/kerem/tasks/longrna/rawdatas/HSC_{sample}.fastq.gz").sample)
    output:
        directory("/mnt/groupMansuy/kerem/tasks/longrna/exp/snakemake/01_FastQC_raw_data")
    conda:
        "long01"
    shell:
        "bash 01_FastQC_raw_data.sh {input} {output}"

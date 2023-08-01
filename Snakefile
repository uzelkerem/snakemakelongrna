configfile: "config.yaml"

rule fastqc_raw_data:
    input:
        "data/{sample}.fastq.gz"
    output:
        "exp/preprocess_01/01_FastQC_raw_data"
    log:
        "logs/01_FastQC/{sample}.log"
    conda:
        "envs/long01.yaml"
    threads: 4
    shell:
        "(fastqc -t {threads} {input} > {output}) 2> {log}"
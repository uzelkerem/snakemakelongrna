configfile: "config.yaml"

rule all:
    input:
        "results/preprocess_01/01_FastQC_raw_data/.dir",
        "results/preprocess_01/02_UMI_extraction/.dir",
        "results/preprocess_01/03_trimmed_data/.dir",
        "results/preprocess_01/04_FastQC_trimmed_data/.dir",
        "results/preprocess_01/05_sortmernaed_data/.dir",
        expand(
            [
                "results/preprocess_01/04_FastQC_trimmed_data/{sample}_R1_processed_trimmed_fastqc.html",
                "results/preprocess_01/05_sortmernaed_data/zipped/{sample}_R1_processed_trimmed_other.fq.gz",
                "results/preprocess_01/01_FastQC_raw_data/{sample}_R1_fastqc.html",
                "results/preprocess_01/01_FastQC_raw_data/{sample}_R2_fastqc.zip"
            ],
            sample=config["samples"]
        )

rule setup_directories:
    output:
        [f"{d}/.dir" for d in config['directories']]
    shell:
        """
        {commands}
        """.format(commands="\n".join(f"mkdir -p {d} && touch {d}/.dir" for d in config['directories']))


rule fqc_raw_data_R1:
    input:
        fastq="data/{sample}_R1.fastq.gz"
    output:
        html="results/preprocess_01/01_FastQC_raw_data/{sample}_R1_fastqc.html",
        zip="results/preprocess_01/01_FastQC_raw_data/{sample}_R1_fastqc.zip"
    log:
        "logs/01_FastQC_rawdata/{sample}_R1.log"
    conda:
        "envs/fastq.yaml"
    threads: 4
    shell:
        "fastqc -t {threads} -o results/preprocess_01/01_FastQC_raw_data/ {input.fastq} > {log} 2>&1"

rule fqc_raw_data_R2:
    input:
        fastq="data/{sample}_R2.fastq.gz"
    output:
        html="results/preprocess_01/01_FastQC_raw_data/{sample}_R2_fastqc.html",
        zip="results/preprocess_01/01_FastQC_raw_data/{sample}_R2_fastqc.zip"
    log:
        "logs/01_FastQC_rawdata/{sample}_R2.log"
    conda:
        "envs/fastq.yaml"
    threads: 4
    shell:
        "fastqc -t {threads} -o results/preprocess_01/01_FastQC_raw_data/ {input.fastq} > {log} 2>&1"

rule umi_extraction:
    wildcard_constraints:
        sample="[^/]+"
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        out_r1="results/preprocess_01/02_UMI_extraction/{sample}_R1_processed.fastq.gz",
        out_r2="results/preprocess_01/02_UMI_extraction/Read2s/{sample}_R2_processed.fastq.gz"
    log:
        "logs/02_UMI_extraction/{sample}_extraction_log.txt"
    conda:
        "envs/umi.yaml"
    shell:
        """
        umi_tools extract -I {input.r2} -S {output.out_r2} \
            --read2-in {input.r1} --read2-out {output.out_r1} \
            --bc-pattern NNNNNNNN \
            --log {log}
        """

rule trimming:
    input:
        fastq="results/preprocess_01/02_UMI_extraction/{sample}_R1_processed.fastq.gz"
    output:
        trimmed="results/preprocess_01/03_trimmed_data/{sample}_R1_processed_trimmed.fq.gz",
        report="results/preprocess_01/03_trimmed_data/{sample}_R1_processed.fastq.gz_trimming_report.txt"
    params:
        quality=config["trimming_params"]["quality"],
        length=config["trimming_params"]["length"],
        stringency=config["trimming_params"]["stringency"]
    log:
        "logs/03_trimming/{sample}_R1_processed.log"
    conda:
        "envs/trimg.yaml"
    shell:
        """
        trim_galore {input.fastq} \
            -q {params.quality} \
            --length {params.length} \
            --trim-n \
            --stringency {params.stringency} \
            -o results/preprocess_01/03_trimmed_data/ \
            > {log} 2>&1
        """

rule fqc_trimmed_data:
    input:
        trimmed_fq="results/preprocess_01/03_trimmed_data/{sample}_R1_processed_trimmed.fq.gz"
    output:
        html="results/preprocess_01/04_FastQC_trimmed_data/{sample}_R1_processed_trimmed_fastqc.html",
        zip="results/preprocess_01/04_FastQC_trimmed_data/{sample}_R1_processed_trimmed_fastqc.zip"
    log:
        "logs/04_FastQC_trimmed/{sample}_R1_trimmed.log"
    conda:
        "envs/fastq.yaml"
    threads: 4
    shell:
        "fastqc -t {threads} -o results/preprocess_01/04_FastQC_trimmed_data/ {input.trimmed_fq} > {log} 2>&1"

rule unzip_trimmed:
    input:
        trimmed_fq="results/preprocess_01/03_trimmed_data/{sample}_R1_processed_trimmed.fq.gz"
    output:
        fastq="results/preprocess_01/03_trimmed_data/unzipped/{sample}_R1_processed_trimmed.fq"
    shell:
        "gunzip -c {input.trimmed_fq} > {output.fastq}"

rule sort_rRNAs_out:
    input:
        trimmed_fq="results/preprocess_01/03_trimmed_data/unzipped/{sample}_R1_processed_trimmed.fq",
        refs=expand("{ref}", ref=config["sortmerna_ref"])
    output:
        aligned="results/preprocess_01/05_sortmernaed_data/{sample}_R1_processed_trimmed_aligned.fq",
        other="results/preprocess_01/05_sortmernaed_data/{sample}_R1_processed_trimmed_other.fq"
    log:
        "logs/05_sortmerna/{sample}_sortmernaed.log"
    conda:
        "envs/sortrna.yaml"
    threads: 8
    shell:
        """
        # Create directories if they don't exist
        working_dir="results/preprocess_01/05_sortmernaed_data/tempdir{wildcards.sample}"
        mkdir -p $working_dir
        
        # Set up the sortmerna command
        ref_str=""
        for ref in {input.refs}; do
            ref_str="$ref_str --ref $ref"
        done
        
        # Run sortmerna
        sortmerna $ref_str \
        --workdir $working_dir \
        --reads {input.trimmed_fq} --threads {threads} \
        --aligned results/preprocess_01/05_sortmernaed_data/{wildcards.sample}_R1_processed_trimmed_aligned \
        --other results/preprocess_01/05_sortmernaed_data/{wildcards.sample}_R1_processed_trimmed_other \
        --fastx > {log} 2>&1
        
        # Optionally clean up the temporary directories if needed
        rm -r $working_dir
        """
        
rule zip_sortmernaed:
    input:
        other_fq="results/preprocess_01/05_sortmernaed_data/{sample}_R1_processed_trimmed_other.fq"
    output:
        fqgz="results/preprocess_01/05_sortmernaed_data/zipped/{sample}_R1_processed_trimmed_other.fq.gz"
    shell:
        "gzip -c {input.other_fq} > {output.fqgz}"  

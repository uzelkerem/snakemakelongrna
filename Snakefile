configfile: "config.yaml"

rule all:
    input:
        "results/preprocess_01/01_FastQC_raw_data/.dir",
        "results/preprocess_01/02_UMI_extraction/.dir",
        "results/preprocess_01/03_trimmed_data/.dir",
        "results/preprocess_01/04_FastQC_trimmed_data/.dir",
        "results/preprocess_01/05_sortmernaed_data/.dir",
        "results/preprocess_01/06_FQScreen/.dir",
        "results/preprocess_01/07_star_aligned/.dir",
        "results/preprocess_01/08_umi_deduplicated/.dir",
        "results/preprocess_01/09_picard_markdup/.dir",
        "results/preprocess_01/10_featureCounts/.dir",
        "results/preprocess_01/11_TinScore/.dir",
        "results/preprocess_01/12_GeneBodyCov/.dir",
        "results/preprocess_01/10_featureCounts/{prefix}_counts_gtfD_s02_sortmerna.txt".format(prefix=config['prefix']),
        "results/preprocess_01/11_TinScore/merged.tsv",
        expand(
            [
                "results/preprocess_01/04_FastQC_trimmed_data/{sample}_R1_processed_trimmed_fastqc.html",
                "results/preprocess_01/06_FQScreen/from_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.html",
                "results/preprocess_01/09_picard_markdup/afterumidedup/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup_marked_duplicates_metrics.txt",
                "results/preprocess_01/12_GeneBodyCov/{sample}.geneBodyCoverage.txt",
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
        fastq=temp("results/preprocess_01/03_trimmed_data/unzipped/{sample}_R1_processed_trimmed.fq")
    shell:
        "gunzip -c {input.trimmed_fq} > {output.fastq}"

rule sort_rRNAs_out:
    input:
        trimmed_fq="results/preprocess_01/03_trimmed_data/unzipped/{sample}_R1_processed_trimmed.fq",
        refs=expand("{ref}", ref=config["sortmerna_ref"])
    output:
        other=temp("results/preprocess_01/05_sortmernaed_data/{sample}_R1_processed_trimmed_other.fq")
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

rule fqscreen_trimmed_data:
    input:
        trimmed_fq="results/preprocess_01/03_trimmed_data/{sample}_R1_processed_trimmed.fq.gz",
        sortmernaed_fq="results/preprocess_01/05_sortmernaed_data/zipped/{sample}_R1_processed_trimmed_other.fq.gz"
    output:
        trimmed_fqscreen_html="results/preprocess_01/06_FQScreen/from_trimmed_data/{sample}_R1_processed_trimmed_screen.html",
        trimmed_fqscreen_txt="results/preprocess_01/06_FQScreen/from_trimmed_data/{sample}_R1_processed_trimmed_screen.txt",
        sortmernaed_fqscreen_html="results/preprocess_01/06_FQScreen/from_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.html",
        sortmernaed_fqscreen_txt="results/preprocess_01/06_FQScreen/from_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.txt"        
    log:
        trimmed="logs/06_FQScreen/{sample}_R1_trimmed.log",
        sortmernaed="logs/06_FQScreen/{sample}_R1_sortmernaed.log"
    conda:
        "envs/fqscreen.yaml"
    threads: 4
    shell:
        """
        fastq_screen --aligner bowtie2 --threads {threads} \
        --conf {config[fastq_screen_conf]} \
        --outdir results/preprocess_01/06_FQScreen/from_trimmed_data {input.trimmed_fq} > {log.trimmed} 2>&1

        fastq_screen --aligner bowtie2 --threads {threads} \
        --conf {config[fastq_screen_conf]} \
        --outdir results/preprocess_01/06_FQScreen/from_sortmernaed_data {input.sortmernaed_fq} > {log.sortmernaed} 2>&1
        """

rule star_aligner:
    input:
        fastq="results/preprocess_01/05_sortmernaed_data/zipped/{sample}_R1_processed_trimmed_other.fq.gz",
    output:
        bam="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.bam",
        idx="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.bam.bai",
        idxstat="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.idxstats",
        flagstat="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.flagstat"
    log:
        star="logs/07_staraligner/{sample}_staralign.log",
        index="logs/07_staraligner/{sample}_index.log",
    conda:
        "envs/star.yaml"
    threads: 8
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {config[star_genome_dir]} \
            --sjdbGTFfile {config[GTF_file]} \
            --readFilesCommand zcat \
            --readFilesIn {input.fastq} \
         --outFileNamePrefix results/preprocess_01/07_star_aligned/{wildcards.sample}_R1_processed_trimmed_other_ \
         --outSAMtype BAM SortedByCoordinate > {log.star} 2>&1

        # Generate index for the output BAM file
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.idxstat} 
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule umi_deduplication:
    input:
        bam="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.bam"
    output:
        dedup_bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam"
    log:
        "logs/08_UMI_deduplication/{sample}_deduplication_log.txt"
    conda:
        "envs/umi.yaml"
    shell:
        """
        # Deduplicate BAM file using umi_tools dedup
        umi_tools dedup \
        --stdin={input.bam} \
        --log={log} \
        --stdout {output.dedup_bam} --buffer-whole-contig

        # Generate index for the deduplicated BAM file
        samtools index {output.dedup_bam}
        """

rule picard_markduplicates:
    input:
        beforeumi_bam="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.bam",
        afterumi_bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam"
    output:
        beforeumi_metrics="results/preprocess_01/09_picard_markdup/beforeumidedup/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out_marked_duplicates_metrics.txt",
        afterumi_metrics="results/preprocess_01/09_picard_markdup/afterumidedup/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup_marked_duplicates_metrics.txt"
    log:
        beforeumi="logs/09_picard_markdup/{sample}_beforeumi_markdup_log.txt",
        afterumi="logs/09_picard_markdup/{sample}_afterumi_markdup_log.txt"
    conda:
        "envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
        I={input.beforeumi_bam} \
        O="/dev/null" \
        M={output.beforeumi_metrics} \
        > {log.beforeumi} 2>&1

        picard MarkDuplicates \
        I={input.afterumi_bam} \
        O="/dev/null" \
        M={output.afterumi_metrics} \
        > {log.afterumi} 2>&1
        """

rule feature_counts:
    input:
        bams=expand("results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam", sample=config["samples"])
    output:
        counts="results/preprocess_01/10_featureCounts/{prefix}_counts_gtfD_s02_sortmerna.txt".format(prefix=config['prefix']),
        summary="results/preprocess_01/10_featureCounts/{prefix}_counts_gtfD_s02_sortmerna.txt.summary".format(prefix=config['prefix'])
    log:
        "logs/10_featureCounts/{prefix}_counts_gtfD_s02_sortmerna.log".format(prefix=config['prefix'])
    conda:
        "envs/featurecounts.yaml"
    threads: 8
    shell:
        """
        featureCounts -T {threads} -s 2 -t exon -g gene_id -a {config[GTF_file]} \
        -o {output.counts} {input.bams} > {log} 2>&1
        """

rule calculate_tin:
    input:
        bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam"
    output:
        tsv="results/preprocess_01/11_TinScore/{sample}.tsv"
    conda:
        "envs/rseqc.yaml"
    threads: 8

    shell:
        "calculate-tin.py -r {config[bed_file]} -i {input.bam} --names={wildcards.sample} -p {threads} 1> {output.tsv}"

rule merge_tin:
    input:
        expand("results/preprocess_01/11_TinScore/{sample}.tsv", sample=config["samples"])
    output:
        "results/preprocess_01/11_TinScore/merged.tsv"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        merge-tin.py --input-files {input} --output-file {output}
        """

rule calculate_genebodycoverage:
    input:
        bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam"
    output:
        pdf="results/preprocess_01/12_GeneBodyCov/{sample}.geneBodyCoverage.curves.pdf",
        r="results/preprocess_01/12_GeneBodyCov/{sample}.geneBodyCoverage.r",
        txt="results/preprocess_01/12_GeneBodyCov/{sample}.geneBodyCoverage.txt"
    conda:
        "envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -r {config[bed_file]} -i {input.bam} -o results/preprocess_01/12_GeneBodyCov/{wildcards.sample}"
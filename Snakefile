configfile: "config.yaml"

# Expand analysis directories and their subdirectories for rule all targets
analysis_targets = []
for analysis in config["ANALYSIS_DIRS"]:
    # Skip if run_rseqc is False and the analysis is 12_GeneBodyCov
    if (not config["run_rseqc"] and analysis == "12_GeneBodyCov"):
        continue
    analysis_targets.append(f"results/qc_plots_02/{analysis}/multiqc_report.html")

rule all:
    input:
        analysis_targets,
        "results/preprocess_01/01_FastQC_raw_data/.dir",
        "results/preprocess_01/02_UMI_extraction/.dir",
        "results/preprocess_01/03_trimmed_data/.dir",
        "results/preprocess_01/04_FastQC_trimmed_data/.dir",
        "results/preprocess_01/05_sortmernaed_data/.dir",
        "results/preprocess_01/06_01_FQScreen_trimmed_data/.dir",
        "results/preprocess_01/06_02_FQScreen_sortmernaed_data/.dir",
        "results/preprocess_01/07_star_aligned/.dir",
        "results/preprocess_01/08_umi_deduplicated/.dir",
        "results/preprocess_01/09_01_markdup_beforeumidedup/.dir",
        "results/preprocess_01/09_02_markdup_afterumidedup/.dir",
        "results/preprocess_01/10_featureCounts/.dir",
        "results/preprocess_01/11_TinScore/.dir",
        "results/preprocess_01/12_GeneBodyCov/.dir",
        "results/qc_plots_02/01_FastQC_raw_data/.dir",
        "results/qc_plots_02/03_trimmed_data/.dir",
        "results/qc_plots_02/04_FastQC_trimmed_data/.dir",
        "results/qc_plots_02/06_01_FQScreen_trimmed_data/.dir",
        "results/qc_plots_02/06_02_FQScreen_sortmernaed_data/.dir",
        "results/qc_plots_02/07_star_aligned/.dir",
        "results/qc_plots_02/08_umi_deduplicated/.dir",
        "results/qc_plots_02/09_01_markdup_beforeumidedup/.dir",
        "results/qc_plots_02/09_02_markdup_afterumidedup/.dir",
        "results/qc_plots_02/10_featureCounts/.dir",
        "results/qc_plots_02/11_TinScore/.dir",
        "results/qc_plots_02/12_GeneBodyCov/.dir",
        "results/preprocess_01/10_featureCounts/{prefix}_counts_gtfD_s02_sortmerna.txt".format(prefix=config['prefix']),
        "results/preprocess_01/11_TinScore/merged.tsv" if config["run_rseqc"] else None,
        expand("results/preprocess_01/12_GeneBodyCov/{sample}.geneBodyCoverage.txt", sample=config["samples"]) if config["run_rseqc"] else [],
        expand(
            [
                "results/preprocess_01/04_FastQC_trimmed_data/{sample}_R1_processed_trimmed_fastqc.html",
                "results/preprocess_01/06_01_FQScreen_trimmed_data/{sample}_R1_processed_trimmed_other_screen.html",
                "results/preprocess_01/06_02_FQScreen_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.html",
                "results/preprocess_01/09_02_markdup_afterumidedup/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup_marked_duplicates_metrics.txt",
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
        out_r2=temp("results/preprocess_01/02_UMI_extraction/Read2s/{sample}_R2_processed.fastq.gz")
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
        report="results/preprocess_01/03_trimmed_data/{sample}_R1_processed_trimming_report.txt"
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
        
        # Rename the trimming report
        mv results/preprocess_01/03_trimmed_data/{wildcards.sample}_R1_processed.fastq.gz_trimming_report.txt {output.report}
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
        trimmed_fqscreen_html="results/preprocess_01/06_01_FQScreen_trimmed_data/{sample}_R1_processed_trimmed_screen.html",
        trimmed_fqscreen_txt="results/preprocess_01/06_01_FQScreen_trimmed_data/{sample}_R1_processed_trimmed_screen.txt",
        sortmernaed_fqscreen_html="results/preprocess_01/06_02_FQScreen_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.html",
        sortmernaed_fqscreen_txt="results/preprocess_01/06_02_FQScreen_sortmernaed_data/{sample}_R1_processed_trimmed_other_screen.txt"        
    log:
        trimmed="logs/06_01_FQScreen_trimmed_data/{sample}_R1_trimmed.log",
        sortmernaed="logs/06_02_FQScreen_sortmernaed_data/{sample}_R1_sortmernaed.log"
    conda:
        "envs/fqscreen.yaml"
    threads: 4
    shell:
        """
        fastq_screen --aligner bowtie2 --threads {threads} \
        --conf {config[fastq_screen_conf]} \
        --outdir results/preprocess_01/06_01_FQScreen_trimmed_data {input.trimmed_fq} > {log.trimmed} 2>&1

        fastq_screen --aligner bowtie2 --threads {threads} \
        --conf {config[fastq_screen_conf]} \
        --outdir results/preprocess_01/06_02_FQScreen_sortmernaed_data {input.sortmernaed_fq} > {log.sortmernaed} 2>&1
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
        dedup_bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam",
        copied_log="results/preprocess_01/08_umi_deduplicated/{sample}_deduplication_log.txt"
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

        # Copy the log file to the desired output directory
        cp {log} {output.copied_log}
        """

rule picard_markduplicates:
    input:
        beforeumi_bam="results/preprocess_01/07_star_aligned/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out.bam",
        afterumi_bam="results/preprocess_01/08_umi_deduplicated/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup.bam"
    output:
        beforeumi_metrics="results/preprocess_01/09_01_markdup_beforeumidedup/{sample}_R1_processed_trimmed_other_Aligned.sortedByCoord.out_marked_duplicates_metrics.txt",
        afterumi_metrics="results/preprocess_01/09_02_markdup_afterumidedup/{sample}_R1_processed_trimmed_other_Aligned_sorted_dedup_marked_duplicates_metrics.txt"
    log:
        beforeumi="logs/09_01_markdup_beforeumidedup/{sample}_beforeumi_markdup_log.txt",
        afterumi="logs/09_02_markdup_afterumidedup/{sample}_afterumi_markdup_log.txt"
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
        """
        if [ "{config[run_rseqc]}" = "True" ]; then
            calculate-tin.py -r {config[bed_file]} -i {input.bam} --names={wildcards.sample} -p {threads} 1> {output.tsv}
        else
            echo "RSeQC analysis is turned off in the config."
            exit 1
        fi
        """

rule merge_tin:
    input:
        expand("results/preprocess_01/11_TinScore/{sample}.tsv", sample=config["samples"])
    output:
        "results/preprocess_01/11_TinScore/merged.tsv"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        if [ "{config[run_rseqc]}" = "True" ]; then
            merge-tin.py --input-files {input} --output-file {output}
        else
            echo "RSeQC analysis is turned off in the config."
            exit 1
        fi
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
        """
        if [ "{config[run_rseqc]}" = "True" ]; then
            geneBody_coverage.py -r {config[bed_file]} -i {input.bam} -o results/preprocess_01/12_GeneBodyCov/{wildcards.sample}
        else
            echo "RSeQC analysis is turned off in the config."
            exit 1
        fi
        """

rule multiqc:
    input:
        lambda wildcards: expand("results/preprocess_01/{analysis}", analysis=wildcards.analysis)
    output:
        html="results/qc_plots_02/{analysis}/multiqc_report.html"
    params:
        outdir="results/qc_plots_02/{analysis}",
        configfile="envs/multiqc_config.yaml"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc -c {params.configfile} -o {params.outdir} {input} -p -s -d
        """

rule plot_tin_scores:
    input:
        "results/preprocess_01/11_TinScore/merged.tsv"
    output:
        "results/qc_plots_02/11_TinScore/tin_scores.png"
    conda:
        "envs/r_ggplot2.yaml"
    shell:
        """
        if [ "{config[run_rseqc]}" = "True" ]; then
            Rscript scripts/plot_tin_scores.R {input} {output}
        else
            echo "RSeQC analysis is turned off in the config."
            exit 1
        fi
        """
#test_glob_wildcards.py

from snakemake import glob_wildcards

input_folder = "/mnt/groupMansuy/kerem/tasks/longrna/rawdatas"

wildcards = glob_wildcards(f"{input_folder}/{{sample}}.fq.gz")
print(wildcards)

sample_values = wildcards.sample
print(sample_values)

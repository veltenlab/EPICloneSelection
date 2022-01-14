configfile: "config_split.yaml"

import pandas as pd
import subprocess

samples = config['samples']

rule all:
    input:
#        expand("../data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",sample=samples),
#        expand("../data/external/Cabezas/raw/{sample}/trimmed/{sample}_1_val_1.fq",sample=samples),
        expand("../data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam",sample=samples)

#rule download:
#    output:
#        folder = "../data/external/Cabezas/raw/{sample}/fastq/",
#        first = "../data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",
#        second = "../data/external/Cabezas/raw/{sample}/fastq/{sample}_2.fastq"
#    params:
#        s = "{sample}",
#        max_size = "50G"
#    log:
#        "../data/external/Cabezas/raw/{sample}/fastq/{sample}.log"
#    run:
#        print("Start " + params.s)
#        sra_download = "/users/lvelten/mscherer/conda/envs/sra/bin/prefetch " + params.s + " --output-directory " + output.folder + " --max-size " + params.max_size
#        proc = subprocess.run(sra_download,shell=True)
#        dump_command = "/users/lvelten/mscherer/conda/envs/sra/bin/fastq-dump --split-3 " + output.folder + "/" + params.s + "/" + params.s + ".sra"+" -O " + output.folder
#        proc = subprocess.run(dump_command,shell=True)

#rule trim:
#    input:
#        first = "../data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",
#        second = "../data/external/Cabezas/raw/{sample}/fastq/{sample}_2.fastq"
#    output:
#        folder = "../data/external/Cabezas/raw/{sample}/trimmed/",
#        first = "../data/external/Cabezas/raw/{sample}/trimmed/{sample}_1_val_1.fq",
#        second = "../data/external/Cabezas/raw/{sample}/trimmed/{sample}_2_val_2.fq"
#    conda:
#        "/users/lvelten/mscherer/conda/envs/genotools.yml"
#    log:
#        "../data/external/Cabezas/raw/{sample}/trimmed/{sample}.log"
#    shell:
#        "trim_galore --nextera --paired {input.first} {input.second} -o {output} ; "
#        "rm -rf {input.first}; rm -rf {input.second} ; "
#        "2> {log}"

rule bismark:
    input:
        ref = "../references/mm10/",
        first = "../data/external/Cabezas/raw/{sample}/trimmed/{sample}_1_val_1.fq",
        second = "../data/external/Cabezas/raw/{sample}/trimmed/{sample}_2_val_2.fq"
    output:
        "../data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam"
    conda:
        "/users/lvelten/mscherer/conda/envs/bismark.yml"
    threads:
        4
    log:
        "../data/external/Cabezas/raw/{sample}/mapped/{sample}.log"
    params:
        threads = 4,
        tmp_dir = "tmp"
    shell:
        "bismark --bowtie2 {input.ref} -1 {input.first} -2 {input.second} -o {output} --temp_dir {params.tmp_dir} -p {params.threads} ; "
        "rm -rf {input.first}; rm -rf {input.second} ; "
        "2> {log}"
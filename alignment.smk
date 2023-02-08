import os
import pandas as pd


THREADS=24
TEMP="/tmp"

DATASET=config["dataset"]
SAMPLES=list(config["samples"].keys())
print("DATASET: %s" % DATASET)
for sample in SAMPLES:
    print("SAMPLE: %s" % sample)
    print("R1: %s" % config["samples"][sample]["R1"])
    print("R2: %s" % config["samples"][sample]["R2"])
R1s=lambda wildcards: config["samples"][wildcards.sample]["R1"],
R2s=lambda wildcards: config["samples"][wildcards.sample]["R2"]

REF=os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta")

rule all:
    input:
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.sa"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.amb"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.ann"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.alt"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.bwt"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.pac"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.dict"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.fai"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta"),

        expand(os.path.join(DATASET, "Stats", "{sample}_fastp.html"), sample=SAMPLES),
        expand(os.path.join(DATASET, "Stats", "{sample}_fastp.json"), sample=SAMPLES),
        expand(os.path.join(DATASET, "Fastp_out", "{sample}_1.fastq"), sample=SAMPLES),
        expand(os.path.join(DATASET, "Fastp_out", "{sample}_2.fastq"), sample=SAMPLES),

        expand(os.path.join(DATASET, "Alignment", "{sample}.cram"), sample=SAMPLES),
        expand(os.path.join(DATASET, "Alignment", "{sample}.cram.crai"), sample=SAMPLES),

        expand(os.path.join(DATASET, "Stats", "{sample}.cram.stats"), sample=SAMPLES),

        expand(os.path.join(DATASET, "Stats", "{sample}.PEFragmentLength.txt"), sample=SAMPLES),
        expand(os.path.join(DATASET, "Stats", "{sample}.PEFragmentLength.png"), sample=SAMPLES),
        expand(os.path.join(DATASET, "PEFragmentSize", "{sample}.PEFragmentLength.tsv"), sample=SAMPLES)


# Get GATK BWA index
if config["reference"] == "":
    rule get_index:
        output:
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.sa"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.amb"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.ann"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.alt"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.bwt"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.pac"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.dict"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.fai"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta"),
        conda:
            "alignment.yaml"
        shell:
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.sa {output[0]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.amb {output[1]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.ann {output[2]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.alt {output[3]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.bwt {output[4]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.pac {output[5]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict {output[6]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai {output[7]} && "
            "aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta {output[8]} "
else:
    rule get_index:
        output:
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.sa"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.amb"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.ann"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.alt"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.bwt"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.pac"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.dict"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.fai"),
            os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta"),
            os.path.join(DATASET, "BWAIndex")
        resources:
            ref_dict=config["reference"]
        shell:
            "mkdir -p {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.sa {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.amb {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.ann {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.alt {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.bwt {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.64.pac {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.dict {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta.fai {output[9]} && "
            "ln -s {resources.ref_dict}/Homo_sapiens_assembly38.fasta {output[9]} "
            

# Get fastq stats with fastp
rule fastp:
    input:
        R1s,
        R2s, 

        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.sa"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.amb"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.ann"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.alt"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.bwt"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.64.pac"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.dict"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta.fai"),
        os.path.join(DATASET, "BWAIndex", "Homo_sapiens_assembly38.fasta")
        # lambda wildcards: config["samples"][wildcards.sample]["R1"],
        # lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        temp(os.path.join(DATASET, "Fastp_out", "{sample}_1.fastq")),
        temp(os.path.join(DATASET, "Fastp_out", "{sample}_2.fastq")),
        os.path.join(DATASET, "Stats", "{sample}_fastp.html"),
        os.path.join(DATASET, "Stats", "{sample}_fastp.json")
    threads:
        16
    resources:
        tmpdir=TEMP
    params:
        DATASET
    conda:
        "alignment.yaml"
    shell:
        "fastp -Q -R \"{params} - {wildcards.sample} Fastp Analysis\" "
        "-i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        "-w {threads} -h {output[2]} -j {output[3]}"


# Run bwa mem in PE mode
rule alignment:
    input:
        os.path.join(DATASET, "Fastp_out", "{sample}_1.fastq"),
        os.path.join(DATASET, "Fastp_out", "{sample}_2.fastq")
    output:
        os.path.join(DATASET, "Alignment", "{sample}.cram")
    threads:
        THREADS
    resources:
        tmpdir=TEMP,
        ref=REF
    conda:
        "alignment.yaml"
    shell:
        "bwa mem -t {threads} -R "
        "\"@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}\\tLB:Seq\" "
        "{resources.ref} {input[0]} {input[1]} | "
        "samblaster -M --ignoreUnmated | "
        "samtools view -@ {threads} -Shb - | "
        "samtools sort -m 1G -O bam -l 0 -T {resources.tmpdir} - | "
        "samtools view -T {resources.ref} -C -o {output} - "


rule index:
    input:
        os.path.join(DATASET, "Alignment", "{sample}.cram")
    output:
        os.path.join(DATASET, "Alignment", "{sample}.cram.crai")
    threads:
        THREADS
    resources:
        tmpdir=TEMP,
        ref=REF
    conda:
        "alignment.yaml"
    shell:
        "samtools index -@ {threads} {input}"


rule cram_stats:
    input:
        os.path.join(DATASET, "Alignment", "{sample}.cram"),
        os.path.join(DATASET, "Alignment", "{sample}.cram.crai")
    output:
        os.path.join(DATASET, "Stats", "{sample}.cram.stats")
    threads:
        THREADS
    resources:
        tmpdir=TEMP,
        ref=REF
    conda:
        "alignment.yaml"
    shell:
        "samtools stats -@ {threads} --reference {resources.ref} {input[0]} > {output}"


rule PEFragmentSize:
    input:
        os.path.join(DATASET, "Alignment", "{sample}.cram"),
        os.path.join(DATASET, "Alignment", "{sample}.cram.crai")
    output:
        os.path.join(DATASET, "Stats", "{sample}.PEFragmentLength.txt"),
        os.path.join(DATASET, "Stats", "{sample}.PEFragmentLength.png"),
        os.path.join(DATASET, "PEFragmentSize", "{sample}.PEFragmentLength.tsv")
    threads:
        THREADS
    resources:
        tmpdir=TEMP,
        ref=REF
    conda:
        "alignment.yaml"
    shell:
        "bamPEFragmentSize -b {input[0]} --table {output[0]} -hist {output[1]} "
        "--maxFragmentLength 1000 --samplesLabel {wildcards.sample} -p {threads} "
        "--outRawFragmentLengths {output[2]} -n 10000 "




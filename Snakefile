
#rule all:
#    input:
#        "fragments/H3K27ac_H3K27ac.tsv.gz",
#        "fragments/H3K27ac_H3K27me3.tsv.gz",
#        "fragments/H3K27ac_RNAPII.tsv.gz",
#        "fragments/H3K27me3_H3K27me3.tsv.gz",
#        "fragments/H3K27me3_H3K27ac.tsv.gz",
#        "fragments/H3K27me3_RNAPII.tsv.gz",
#        "fragments/RNAPII_RNAPII.tsv.gz",
#        "fragments/RNAPII_H3K27ac.tsv.gz",
#        "fragments/RNAPII_H3K27me3.tsv.gz"

rule get_genome:
    output: "genome/mm10.fa.gz"
    threads: 1
    message: "Download mm10 genome"
    shell:
        """
        cd genome
        wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
        """

rule bwa_build:
    input: "genome/mm10.fa.gz"
    output: "genome/mm10.fa"
    threads: 1
    message: "Build bwa index for genome"
    shell:
        """
        gzip -d genome/mm10.fa.gz
        bwa-mem2 index genome/mm10.fa
        """

rule download:
    input: "data/downloads.txt"
    output: "data/SC_AM_index.fastq.gz",
            "data/SC_AM_R1.fastq.gz",
            "data/SC_AM_R2.fastq.gz"
    shell:
        """
        while read line; do
            aws s3 cp $line ./data --request-payer
        done < {input}
        """

rule decompress_fastq:
    input: "data/SC_AM_index.fastq.gz",
           "data/SC_AM_R1.fastq.gz",
           "data/SC_AM_R2.fastq.gz"
    output: "data/SC_AM_index.fastq",
           "data/SC_AM_R1.fastq",
           "data/SC_AM_R2.fastq"
    shell:
        """
        gzip -d {input}
        """

rule attach_barcodes:
    input:
        index="data/SC_AM_index.fastq",
        r1="data/SC_AM_R1.fastq",
        r2="data/SC_AM_R2.fastq"
    output: "data/SC_AM_R1.barcoded.fastq", "data/SC_AM_R2.barcoded.fastq"
    shell:
        """
        sinto barcode \
                --barcode_fastq {input.index} \
                --read1 {input.r1} \
                --read2 {input.r2} \
                --bases 16
        """

#rule split:
#    input: 
#        bc1="data/barcodes_a.fa",
#        bc2="data/barcodes_b.fa",
#    output:
#    threads: 1
#    shell:
#        """
#        # split by tn5 barcode, trim reads
#        """
#
#rule map:
#    input:
#        genome="genome/mm10.fa",
#        r1=,
#        r2=
#    output:
#    threads: 12
#    shell:
#        """
#        bwa-mem2 mem
#        """
#
#rule sort_bam:
#    input:
#    output:
#    threads: 6
#    shell:
#        """
#        samtools sort -@ {threads} -O bam {input} > {output}
#        """
#
#rule fragment:
#    input:
#    output:
#    threads: 12
#    shell:
#        """
#        # sinto fragments
#        """

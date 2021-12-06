rule all:
    input: directory("fragments")

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
    message: "Download raw data from AWS"
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
    message: "Decompress fastq files"
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
    message: "Attach cell barcodes"
    shell:
        """
        sinto barcode \
                --barcode_fastq {input.index} \
                --read1 {input.r1} \
                --read2 {input.r2} \
                --bases 16
        """

rule split:
    input: 
        bc1="data/barcodes_a.fa",
        bc2="data/barcodes_b.fa",
        r1="data/SC_AM_R1.barcoded.fastq",
        r2="data/SC_AM_R2.barcoded.fastq"
    output: directory("demux")
    message: "Split by Tn5 barcode, trim reads"
    threads: 1
    shell:
        """
        python code/demux.py \
            --read1 {input.r1} \
            --read2 {input.r2} \
            --tn5_i5 {input.bc1} \
            --tn5_i7 {input.bc2} \
            --output split
        """

rule map:
    input:
        dir=directory("demux"),
        genome="genome/mm10.fa"
    output: directory("mapped")
    message: "Map reads to genome"
    threads: 24
    shell:
        """
        mkdir {output}
        cd {input.dir}
        for R1 in $(ls -d *.R1.fastq); do
            fname=${{R1%.R1.fastq}}
            R2=$fname.R2.fastq
            bwa-mem2 mem -t {threads} {input.genome} $R1 $R2 \
                | samtools sort -@ {threads} -O bam - \
                > "../{output}/"$fname".bam"
            samtools index -@ {threads} "../{output}/"$fname".bam"
        done
        """

rule fragments:
    input: directory("mapped")
    output: directory("fragments")
    message: "Create fragment file"
    threads: 12
    shell:
        """
        mkdir {output}
        cd {input}
        for bam in $(ls -d *.bam); do
            fname=${{bam%.bam}}
            sinto fragments -p {threads} -b $bam --barcode_regex "[^:]*" -f ../{output}/$fname.tmp
            sort -k1,1 -k2,2n ../{output}/$fname.tmp > ../{output}/$fname.tsv
            bgzip -@ {threads} ../{output}/$fname.tsv
            tabix -p bed ../{output}/$fname.tsv.gz
            rm ../{output}/$fname.tmp
        done
        """
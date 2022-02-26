import pandas as pd
m = pd.read_csv("samples/ronca_metadata_final.csv", header = 0)
PLATES = m['plate_barcode'].unique().tolist() 
SAMPLES = m['well_barcodefull'].unique().tolist() # well barcode
LANES = m['seqsomm'].unique().tolist() # somm
READS = ["1", "2"]
TMPDIR = "/scratch/rapeek"

#print(PLATES)
#print(SAMPLES)
#print(LANES)

rule all:
    input: 
        expand("outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam", lane = LANES, plate = PLATES, sample = SAMPLES)
# should look like this: SOMM504_CGTCTT_R1_GGGCTAACGATGCAGG.fastq after perl output

# remove expand here so that it runs rule once instead twice (for each R1 and R2)
rule unzip:
    input: "../../ronca/raw/{lane}_CKDL200163818-1a_HCJKCCCX2_L7_{read}.fq.gz"
    output: "inputs/fastq/{lane}_R{read}.fastq"
    threads: 1
    resources:
        mem_mb=2000,
	#tmpdir=TMPDIR,
        time=2880
    benchmark: "benchmarks/unzip_fastq_{lane}_R{read}.tsv"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule plate_split_fastq:
    input: "inputs/fastq/{lane}_R{read}.fastq"
    output: "outputs/fastq_plate/{lane}_{plate}_R{read}.fastq"
    threads: 1
    resources:
        mem_mb=2000,
	#tmpdir=TMPDIR,
        time=2880
    #benchmark: "benchmarks/plate_split_{lane}_{plate}_R{read}.tsv" 
    shell:"""
    grep --no-group-separator -A 3 ":{wildcards.plate}" {input} > {output}
    """
# could use checkpoint: flexible, allows unknown number of outputs or output names, wrap directory
rule well_split_fastq:
    input: expand("outputs/fastq_plate/{{lane}}_{{plate}}_R{read}.fastq", read = READS)
    output: expand("outputs/fastq_split/{{lane}}_{{plate}}_R{read}_{sample}.fastq", sample = SAMPLES, read = READS)
    threads: 1
    resources:
        mem_mb=2000,
	#tmpdir=TMPDIR,
        time=2880
    #benchmark: "benchmarks/well_split_fastq_{lane}_{plate}_R{read}_{sample}.tsv"
    params: outdir = "outputs/fastq_split/"
    shell:"""
    /home/rapeek/projects/SEQS/code/BarcodeSplit_RAD_PE.2019.pl {input} GGACAAGCTATGCAGG,GGAAACATCGTGCAGG,GGACATTGGCTGCAGG,GGACCACTGTTGCAGG,GGAACGTGATTGCAGG,GGCGCTGATCTGCAGG,GGCAGATCTGTGCAGG,GGATGCCTAATGCAGG,GGAACGAACGTGCAGG,GGAGTACAAGTGCAGG,GGCATCAAGTTGCAGG,GGAGTGGTCATGCAGG,GGAACAACCATGCAGG,GGAACCGAGATGCAGG,GGAACGCTTATGCAGG,GGAAGACGGATGCAGG,GGAAGGTACATGCAGG,GGACACAGAATGCAGG,GGACAGCAGATGCAGG,GGACCTCCAATGCAGG,GGACGCTCGATGCAGG,GGACGTATCATGCAGG,GGACTATGCATGCAGG,GGAGAGTCAATGCAGG,GGAGATCGCATGCAGG,GGAGCAGGAATGCAGG,GGAGTCACTATGCAGG,GGATCCTGTATGCAGG,GGATTGAGGATGCAGG,GGCAACCACATGCAGG,GGCAAGACTATGCAGG,GGCAATGGAATGCAGG,GGCACTTCGATGCAGG,GGCAGCGTTATGCAGG,GGCATACCAATGCAGG,GGCCAGTTCATGCAGG,GGCCGAAGTATGCAGG,GGCCGTGAGATGCAGG,GGCCTCCTGATGCAGG,GGCGAACTTATGCAGG,GGCGACTGGATGCAGG,GGCGCATACATGCAGG,GGCTCAATGATGCAGG,GGCTGAGCCATGCAGG,GGCTGGCATATGCAGG,GGGAATCTGATGCAGG,GGGACTAGTATGCAGG,GGGAGCTGAATGCAGG,GGGATAGACATGCAGG,GGGCCACATATGCAGG,GGGCGAGTAATGCAGG,GGGCTAACGATGCAGG,GGGCTCGGTATGCAGG,GGGGAGAACATGCAGG,GGGGTGCGAATGCAGG,GGGTACGCAATGCAGG,GGGTCGTAGATGCAGG,GGGTCTGTCATGCAGG,GGGTGTTCTATGCAGG,GGTAGGATGATGCAGG,GGTATCAGCATGCAGG,GGTCCGTCTATGCAGG,GGTCTTCACATGCAGG,GGTGAAGAGATGCAGG,GGTGGAACAATGCAGG,GGTGGCTTCATGCAGG,GGTGGTGGTATGCAGG,GGTTCACGCATGCAGG,GGACACGAGATGCAGG,GGAAGAGATCTGCAGG,GGAAGGACACTGCAGG,GGAATCCGTCTGCAGG,GGAATGTTGCTGCAGG,GGACACTGACTGCAGG,GGACAGATTCTGCAGG,GGAGATGTACTGCAGG,GGAGCACCTCTGCAGG,GGAGCCATGCTGCAGG,GGAGGCTAACTGCAGG,GGATAGCGACTGCAGG,GGACGACAAGTGCAGG,GGATTGGCTCTGCAGG,GGCAAGGAGCTGCAGG,GGCACCTTACTGCAGG,GGCCATCCTCTGCAGG,GGCCGACAACTGCAGG,GGAGTCAAGCTGCAGG,GGCCTCTATCTGCAGG,GGCGACACACTGCAGG,GGCGGATTGCTGCAGG,GGCTAAGGTCTGCAGG,GGGAACAGGCTGCAGG,GGGACAGTGCTGCAGG,GGGAGTTAGCTGCAGG,GGGATGAATCTGCAGG,GGGCCAAGACTGCAGG {params.outdir}{wildcards.lane}_{wildcards.plate}
    """

# rule to align and combine
rule align_fastq:
    input: 
        fq = expand("outputs/fastq_split/{{lane}}_{{plate}}_R{read}_{{sample}}.fastq", read = READS),
	ref = "/home/rapeek/projects/SEQS/final_contigs_300.fa"
    output: "outputs/bams/{lane}_{plate}_{sample}.sort.bam"
    conda: "envs/samtools_bwa.yml"
    threads: 1
    resources:
        mem_mb=2000,
	#tmpdir=TMPDIR,
        time=2880
    #benchmark: "benchmarks/align_fastq_{lane}_{plate}_{sample}.tsv"
    shell:"""
        bwa mem {input.ref} {input.fq} | samtools view -Sb - | samtools sort - -o {output}
	"""

rule filter_bams:
    input: "outputs/bams/{lane}_{plate}_{sample}.sort.bam"
    output: "outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam"
    conda: "envs/samtools_bwa.yml"
    threads: 1
    resources:
        mem_mb=2000,
	#tmpdir=TMPDIR,
        time=2880
    #benchmark: "benchmarks/filter_bams_{lane}_{plate}_{sample}.tsv"
    shell:"""
        samtools view -f 0x2 -b {input} | samtools rmdup - {output}
        """

rule index_bams:
    input: "outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam"
    output: "outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam.bai"
    conda: "envs/samtools_bwa.yml"
    threads: 1
    resources:
        mem_mb=2000,
        time=2880
    shell:"""
        samtools index {input}
	"""

rule bam_stats:
	input:
		bam = "outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam",
		bai = "outputs/bams/{lane}_{plate}_{sample}.sort.flt.bam"
	output:	"outputs/stats/{lane}_{plate}_{sample}.sort.flt.bam.stats"
	threads: 1
	conda: "envs/samtools_bwa.yml"
	shell:"""
	    samtools stats {input.bam} | grep ^SN | cut -f 2- > {output}
	    """

# to add:
# pca w angsd (need angsd enviro)
# sfs
# thetas
# admixture
# fst



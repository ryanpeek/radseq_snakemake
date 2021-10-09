import pandas as pd
m = pd.read_csv("samples/ronca_metadata_final.csv", header = 0)
PLATES = m['plate_barcode'].unique().tolist() 
#SAMPLES = m['well_barcode'].unique().tolist() # well barcode
LANES = m['seqsomm'].unique().tolist() # somm
READS = ["1", "2"]

rule all:
    input: 
        expand("outputs/split_fastq/{lane}_{plate}_R{reads}.fastq", lane = LANES, plate = PLATES, read = READS)

rule unzip:
    input: "../../ronca/raw/{lane}_CKDL200163818-1a_HCJKCCCX2_L7_{read}.fq.gz"
    output: "inputs/fastq/{lane}_R{read}.fastq"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule barcode_split_fastq:
    input: "inputs/fastq/{lane}_R{read}.fastq"
    output: "outputs/split_fastq/{lane}_{plate}_R{read}.fastq" 
    shell:"""
    grep --no-group-separator -A 3 ":{wildcards.plate}" {input} > {output}
    """

#rule well_split_fastq:
#    input: "inputs/split_fastq/{lane}_{plate}_R{read}.fastq"
#    output: "outputs/split_fastq/{lane}_{plate}_{sample}_R{read}.fastq" 
#    shell:"""
#    /home/rapeek/projects/SEQS/code/BarcodeSplit_RAD_PE.2019.pl {output}_R{read}.fastq {output}_R{read}.fastq GGACAAGCTATGCAGG,GGAAACATCGTGCAGG,GGACATTGGCTGCAGG,GGACCACTGTTGCAGG,GGAACGTGATTGCAGG,GGCGCTGATCTGCAGG,GGCAGATCTGTGCAGG,GGATGCCTAATGCAGG,GGAACGAACGTGCAGG,GGAGTACAAGTGCAGG,GGCATCAAGTTGCAGG,GGAGTGGTCATGCAGG,GGAACAACCATGCAGG,GGAACCGAGATGCAGG,GGAACGCTTATGCAGG,GGAAGACGGATGCAGG,GGAAGGTACATGCAGG,GGACACAGAATGCAGG,GGACAGCAGATGCAGG,GGACCTCCAATGCAGG,GGACGCTCGATGCAGG,GGACGTATCATGCAGG,GGACTATGCATGCAGG,GGAGAGTCAATGCAGG,GGAGATCGCATGCAGG,GGAGCAGGAATGCAGG,GGAGTCACTATGCAGG,GGATCCTGTATGCAGG,GGATTGAGGATGCAGG,GGCAACCACATGCAGG,GGCAAGACTATGCAGG,GGCAATGGAATGCAGG,GGCACTTCGATGCAGG,GGCAGCGTTATGCAGG,GGCATACCAATGCAGG,GGCCAGTTCATGCAGG,GGCCGAAGTATGCAGG,GGCCGTGAGATGCAGG,GGCCTCCTGATGCAGG,GGCGAACTTATGCAGG,GGCGACTGGATGCAGG,GGCGCATACATGCAGG,GGCTCAATGATGCAGG,GGCTGAGCCATGCAGG,GGCTGGCATATGCAGG,GGGAATCTGATGCAGG,GGGACTAGTATGCAGG,GGGAGCTGAATGCAGG,GGGATAGACATGCAGG,GGGCCACATATGCAGG,GGGCGAGTAATGCAGG,GGGCTAACGATGCAGG,GGGCTCGGTATGCAGG,GGGGAGAACATGCAGG,GGGGTGCGAATGCAGG,GGGTACGCAATGCAGG,GGGTCGTAGATGCAGG,GGGTCTGTCATGCAGG,GGGTGTTCTATGCAGG,GGTAGGATGATGCAGG,GGTATCAGCATGCAGG,GGTCCGTCTATGCAGG,GGTCTTCACATGCAGG,GGTGAAGAGATGCAGG,GGTGGAACAATGCAGG,GGTGGCTTCATGCAGG,GGTGGTGGTATGCAGG,GGTTCACGCATGCAGG,GGACACGAGATGCAGG,GGAAGAGATCTGCAGG,GGAAGGACACTGCAGG,GGAATCCGTCTGCAGG,GGAATGTTGCTGCAGG,GGACACTGACTGCAGG,GGACAGATTCTGCAGG,GGAGATGTACTGCAGG,GGAGCACCTCTGCAGG,GGAGCCATGCTGCAGG,GGAGGCTAACTGCAGG,GGATAGCGACTGCAGG,GGACGACAAGTGCAGG,GGATTGGCTCTGCAGG,GGCAAGGAGCTGCAGG,GGCACCTTACTGCAGG,GGCCATCCTCTGCAGG,GGCCGACAACTGCAGG,GGAGTCAAGCTGCAGG,GGCCTCTATCTGCAGG,GGCGACACACTGCAGG,GGCGGATTGCTGCAGG,GGCTAAGGTCTGCAGG,GGGAACAGGCTGCAGG,GGGACAGTGCTGCAGG,GGGAGTTAGCTGCAGG,GGGATGAATCTGCAGG,GGGCCAAGACTGCAGG {plate}
#    """

# TODO combine across lanes, use this approach as a template to cat either fastqs or bams
#rule cat_libraries_R2:
#    input: expand("inputs/raw/{sample}_2.fastq.gz", sample = SAMPLES)
#    output: expand("inputs/cat/{library}_2.fastq.gz", library = LIBRARIES)
#    threads: 1
#    resources:
#        mem_mb=4000
#    run: 
#        merge_df = m[['library_name','run_accession']]
#        merge_df = copy.deepcopy(merge_df)
#        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_2.fastq.gz")
#        merge_dict = merge_df.groupby('library_name')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
#        for library in merge_dict.keys():
#            # merge SRR files
#            to_merge = merge_dict[library]
#            # Check if the merged file results from a single or multiple fastq files.
#            # For n-to-1 merging, concatenate input files to produce the output file
#            merge_nb = len(to_merge)
#            if merge_nb > 1:
#                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + library + "_2.fastq.gz"
#            else:
#                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + library + "_2.fastq.gz"
#            os.system(cmd)



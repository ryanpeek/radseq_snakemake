import pandas as pd

m = pd.read_csv("samples/pulex_samples.txt", sep = "\t", header = 0)
SAMPLES = m['Sample'].unique().tolist()
LANE = m['Lane'].unique().tolist() # sommX
READS = ["1", "2"]

rule all:
    input: 
        expand("outputs/split_fastq/{lane}_{sample}_{reads}.fastq", lane = LANE, sample = SAMPLES, reads = READS)


rule unzip:
    input: "inputs/raw/{lane}_{reads}.fastq.gz"
    output: "inputs/fastq/{lane}_{reads}.fastq"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule barcode_split_fastq:
    input: "inputs/fastq/{lane}_{reads}.fastq"
    output: "outputs/split_fastq/{lane}_{sample}_{reads}.fastq" 
    shell:"""
    grep --no-group-separator -A 3 ":{wildcards.sample}" {input} > {output}
    """

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



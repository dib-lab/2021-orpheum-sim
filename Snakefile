import pandas as pd

refseq = pd.read_csv("inputs/refseq_not_in_gtdb_metadata.tsv", sep = "\t", header = 0)
ACC = metadata['assembly_accession'].to_list()

ORPHEUM_DB = ["gtdb-rs202"]
# set constrained k sizes
#dayhoff_ksizes = [14, 16, 18]
protein_ksizes = [10]


rule download_assemblies:
    output: "inputs/assemblies/{acc}_genomic.fna.gz",
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    run:
        row = metadata.loc[metadata['assembly_accession'] == wildcards.acc]
        assembly_ftp = row['ftp_path'].values
        assembly_ftp = assembly_ftp[0]
        assembly_ftp = assembly_ftp + "/*genomic.fna.gz"
        shell("wget -O {output} {assembly_ftp}")

rule gunzip_assemblies:
    input: "inputs/assemblies/{acc}_genomic.fna.gz",
    output: "inputs/assemblies/{acc}_genomic.fna",
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    shell:'''
    gunzip -c {input} > {output}
    '''

rule download_gff:
    output: "inputs/assemblies/{acc}_genomic.gff.gz",
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    run:
        row = metadata.loc[metadata['assembly_accession'] == wildcards.acc]
        assembly_ftp = row['ftp_path'].values
        assembly_ftp = assembly_ftp[0]
        assembly_ftp = assembly_ftp + "/*genomic.gff.gz"
        shell("wget -O {output} {assembly_ftp}")

rule filter_gff_to_cds:
    input: gff="inputs/assemblies/{acc}_genomic.gff.gz",
    output: cds="outputs/cds_gff/{acc}_cds.gff"
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    conda: "envs/rtracklayer.yml"
    script: "scripts/filter_gff_to_cds.R"

rule extract_cds_from_genome:
    input:
        fna="inputs/assemblies/{acc}_genomic.fna",
        cds="outputs/cds_gff/{acc}_cds.gff"
    output: "outputs/cds_fasta/{acc}_cds.fasta"
    threads: 1
    resources: 
        mem_mb=4000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    shell:'''
    bedtools getfasta -fi {input.fna} -bed {input.cds} -name -s > {output}
    '''

rule define_genome_regions_as_bed:
    input: "inputs/assemblies/{acc}_genomic.fna",
    output: "outputs/genome_regions_bed/{acc}_region.sizes"
    threads: 1
    resources: 
        mem_mb=4000,
        tmpdir=TMPDIR
    conda: "envs/samtools.yml"
    shell:'''    
    samtools faidx {input}
    cut -f 1,2 {input}.fai > {output}
    '''

rule sort_cds_gff:
    input:"outputs/cds_gff/{acc}_cds.gff"
    output:"outputs/cds_gff/{acc}_cds_sorted.gff"
    threads: 1
    resources: 
        mem_mb=4000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    shell:'''
    sortBed -i {input} > {output}
    '''

bedtools complement -i GCF_020520145.1_filtered_sorted.gff -g region.sizes > GCF_020520145.1_noncds.bed
bedtools getfasta -fi GCF_020520145.1_ASM2052014v1_genomic.fna -bed GCF_020520145.1_noncds.bed -name -s > GCF_020520145.1_noncds.fna

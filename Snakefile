import pandas as pd

metadata = pd.read_csv("inputs/all_genomes_genbank_info_metadata.tsv", sep = "\t", header = 0)
ACC = metadata['assembly_accession'].to_list()
SEQ = ['cds', 'noncds']

ORPHEUM_DB = ["gtdb-rs202"]
#dayhoff_ksizes = [14, 16, 18]
protein_ksizes = [10]
ALPHA_KSIZE = expand('protein-k{k}', k=protein_ksizes)


TMPDIR = "/scratch/tereiter/"

rule all:
    input:
        expand("outputs/orpheum/{orpheum_db}/{alpha_ksize}/{acc}_{seq}.summary.json", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, acc = ACC, seq = SEQ)

rule download_assemblies:
    output: "inputs/assemblies/{acc}_genomic.fna.gz",
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    run:
        row = metadata.loc[metadata['assembly_accession'] == wildcards.acc]
        assembly_ftp = row['ftp_path_full_genome'].values
        assembly_ftp = assembly_ftp[0]
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
        gff_ftp = row['ftp_path_full_gff'].values
        gff_ftp = gff_ftp[0]
        shell("wget -O {output} {gff_ftp}")

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

rule complement_cds_gff:
    input:
        gff="outputs/cds_gff/{acc}_cds_sorted.gff",
        sizes= "outputs/genome_regions_bed/{acc}_region.sizes"
    output:"outputs/noncds_bed/{acc}_noncds.bed"
    threads: 1
    resources: 
        mem_mb=4000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    shell:'''
    bedtools complement -i {input.gff} -g {input.sizes} > {output}
    '''

rule extract_noncds_from_genome:
    input:
        fna="inputs/assemblies/{acc}_genomic.fna",
        bed="outputs/noncds_bed/{acc}_noncds.bed"
    output:"outputs/noncds_fasta/{acc}_noncds.fasta"
    threads: 1
    resources: 
        mem_mb=4000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    shell:'''
    bedtools getfasta -fi {input.fna} -bed {input.bed} -name -s > {output}
    '''

rule polyester_simulate_reads:
    input: fasta = "outputs/{seq}_fasta/{acc}_{seq}.fasta"
    output: "outputs/polyester/{acc}_{seq}/sample_01.fasta.gz"
    params:
        outdir = lambda wildcards: "outputs/polyester/" + wildcards.acc + "_" + wildcards.seq + "/",
        num_reps = 1,
        read_length = 150,
        simulate_paired = False,
        num_reads_per_transcript=100,
    benchmark: "benchmarks/polyester-{acc}-{seq}.txt"
    threads: 1
    resources:
      mem_mb=16000,
      tmpdir=TMPDIR
    conda: "envs/polyester.yml"
    script: "scripts/polyester.R"

# orpheum index cp'd over from @bluegenes 2021-rank-compare
rule orpheum_translate_reads:
    input: 
        ref="inputs/orpheum_index/{orpheum_db}.{alphabet}-k{ksize}.nodegraph",
        fasta="outputs/polyester/{acc}_{seq}/sample_01.fasta.gz"
    output:
        pep="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_{seq}.coding.faa",
        nuc="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_{seq}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_{seq}.nuc_noncoding.fna",
        csv="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_{seq}.coding_scores.csv",
        json="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_{seq}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum-translate-{acc}-{seq}-{orpheum_db}-{alphabet}-k{ksize}.txt"
    resources:  mem_mb=500000
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

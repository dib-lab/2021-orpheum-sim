import pandas as pd

metadata = pd.read_csv("inputs/all_genomes_genbank_info_metadata.tsv", sep = "\t", header = 0)
ACC = metadata['assembly_accession'].to_list()
SEQ = ['cds', 'noncds']

ORPHEUM_DB = ["gtdb-rs202", "s__Thermosipho_affectus", "g__Thermosipho",
              "f__Fervidobacteriaceae", "o__Thermotogales", "c__Thermotogae", 
              "p__Thermotogota"]

ACC_W_GFF = ['GCA_000008025.1', 'GCA_000012125.1', 'GCA_000012285.1', 'GCA_000020225.1',
             'GCA_000145825.2', 'GCA_000154205.1', 'GCA_000167435.2', 'GCA_000196555.1',
             'GCA_000299235.1', 'GCA_000830885.1', 'GCA_000970205.1', 'GCA_001700755.2',
             'GCA_001884725.2', 'GCA_002006445.1', 'GCA_002442595.2', 'GCA_004006635.1',
             'GCA_006742205.1', 'GCA_009156025.2', 'GCA_013456555.2', 'GCA_015356815.2',
             'GCA_018205295.1', 'GCA_018282115.1', 'GCA_018336995.1', 'GCA_018398935.1', 
             'GCA_018588215.1', 'GCA_019173545.1', 'GCA_019175305.1', 'GCA_019599295.1',
             'GCA_019688735.1', 'GCA_020520145.1', 'GCA_900083515.1', 'GCA_900156205.1',
             'GCA_900478295.1']

#dayhoff_ksizes = [14, 16, 18]
protein_ksizes = [10]
ALPHA_KSIZE = expand('protein-k{k}', k=protein_ksizes)


TMPDIR = "/scratch/tereiter/"

rule all:
    input:
        "outputs/gtdbtk/gtdbtk.bac120.summary.tsv",
        expand("inputs/gffs/{acc_w_gff}_genomic.gff.gz", acc_w_gff = ACC_W_GFF),
        expand("outputs/cds_fasta_paladin/{acc}_cds.sam", acc = ACC),
        expand("outputs/orpheum/{orpheum_db}/{alpha_ksize}/{acc}_{seq}.summary.json", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, acc = ACC, seq = SEQ),
        expand("outputs/orpheum/{orpheum_db}/{alpha_ksize}/{acc}_cds.nuc_noncoding.cut.dedup.only.fna.gz", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, acc = ACC),
        expand("outputs/orpheum_cutoffs/{orpheum_db}/{alpha_ksize}/{acc}.tsv", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, acc = ACC)

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

rule download_gtdbtk_db:
    output: "inputs/gtdbtk_data.tar.gz"
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    shell:'''
    wget -O {output} https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
    '''

rule decompress_gtdbtk_db:
    input:"inputs/gtdbtk_data.tar.gz"
    output: "inputs/release202/manifest.tsv"
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    shell:'''
    tar xf {input} -C inputs/
    '''

rule gtdb_classify_assemblies:
    input:
        fna=expand("inputs/assemblies/{acc}_genomic.fna", acc = ACC),
        gtdb="inputs/release202/manifest.tsv"
    output: "outputs/gtdbtk/gtdbtk.bac120.summary.tsv"
    threads: 8
    resources: 
        mem_mb=264000,
        tmpdir=TMPDIR
    params:
        indir="inputs/assemblies",
        outdir="outputs/gtdbtk",
        dbdir="inputs/release202",
        ext="fna"
    conda: "envs/gtdbtk.yml"
    shell:'''
    real=$(realpath {params.dbdir})
    export GTDBTK_DATA_PATH=${{real}}
    gtdbtk classify_wf --cpus {threads} --genome_dir {params.indir} --out_dir {params.outdir} -x {params.ext}
    '''

# many of the assemblies don't have gff files; standardize by annotating with bakta
# but use these files to assess pseudo genes via NCBI's PGAP pseudo annotation.
rule download_gff:
    output: "inputs/gffs/{acc_w_gff}_genomic.gff.gz",
    threads: 1
    resources: 
        mem_mb=1000,
        tmpdir=TMPDIR
    run:
        row = metadata.loc[metadata['assembly_accession'] == wildcards.acc]
        gff_ftp = row['ftp_path_full_gff'].values
        gff_ftp = gff_ftp[0]
        shell("wget -O {output} {gff_ftp}")

# uncomment to download db; use db that was downloaded for another project
# rule download_bakta_db:
#    output: "inputs/bakta_db/db/version.json"
#    threads: 1
#    resources: mem_mb = 4000
#    conda: "envs/bakta.yml"
#    shell:'''
#    bakta_db download --output {output}
#    '''

rule bakta_assemblies:
    input: 
        fna= "inputs/assemblies/{acc}_genomic.fna",
        db="/home/tereiter/github/2021-orpheum-refseq/inputs/bakta_db/db/version.json"
    output: 
        "outputs/bakta_assemblies/{acc}.faa",
        "outputs/bakta_assemblies/{acc}.gff3",
        "outputs/bakta_assemblies/{acc}.fna",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    benchmark: "benchmarks/bakta_{acc}.txt"
    conda: 'envs/bakta.yml'
    params: 
        dbdir="~/github/2021-orpheum-refseq/inputs/bakta_db/db/",
        outdir = 'outputs/bakta_assemblies/',
    threads: 1
    shell:'''
    bakta --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdir} \
        --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''


rule filter_gff_to_cds:
    input: gff="outputs/bakta_assemblies/{acc}.gff3",
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
    cut -f 1,2 {input}.fai | sort -n -k1 > {output}
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
    benchmark: "benchmarks/polyester-{acc}-{seq}.txt"
    threads: 1
    resources:
      mem_mb=16000,
      tmpdir=TMPDIR
    conda: "envs/polyester.yml"
    script: "scripts/polyester.R"

rule translate_cds_fasta:
    input: "outputs/cds_fasta/{acc}_cds.fasta"
    output: "outputs/cds_fasta_aa/{acc}_cds.faa"
    threads: 1
    resources:
      mem_mb=4000,
      tmpdir=TMPDIR
    conda: "envs/emboss.yml"
    shell:'''
    transeq {input} {output}
    '''
    
rule paladin_index_cds_fasta:
    input: "outputs/cds_fasta_aa/{acc}_cds.faa"
    output: "outputs/cds_fasta_aa/{acc}_cds.faa.pro"
    threads: 1
    resources:
      mem_mb=4000,
      tmpdir=TMPDIR
    conda: "envs/paladin.yml"
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align_polyester_cds_to_determine_orf:
    input: 
        ref="outputs/cds_fasta_aa/{acc}_cds.faa",
        idx="outputs/cds_fasta_aa/{acc}_cds.faa.pro",
        reads="outputs/polyester/{acc}_cds/sample_01.fasta.gz"
    output: "outputs/cds_fasta_paladin/{acc}_cds.sam"
    threads: 1
    resources:
      mem_mb=4000,
      tmpdir=TMPDIR
    conda: "envs/paladin.yml"
    shell:'''
    paladin align -C -t 1 {input.ref} {input.reads} > {output}
    '''

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
    resources:  
        mem_mb=10000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

rule cut_dedup_nuc_noncoding_read_names_cds_only:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.nuc_noncoding.fna",
    output:  "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.nuc_noncoding.cut.dedup.fna.gz",
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    sed '/^>/ s/__.*//g' {input} | awk '/^>/{{f=!d[$1];d[$1]=1}}f' | gzip > {output}
    '''

rule grab_cut_dedup_coding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.coding.faa", 
    output: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.aa_names.cut.dedup.txt",  
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    grep ">" {input} | sed '/^>/ s/__.*//g' | awk '/^>/{{f=!d[$1];d[$1]=1}}f' > {output}
    '''

rule isolate_noncoding_only_reads:
    input: 
        pep ="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.aa_names.cut.dedup.txt",  
        noncoding = "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.nuc_noncoding.cut.dedup.fna.gz",
    output: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{acc}_cds.nuc_noncoding.cut.dedup.only.fna.gz",
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 1
    shell:"""
    zcat {input.noncoding} | paste - - | grep -v -F -f {input.pep} | tr "\t" "\n" | gzip > {output}
    """

rule calculate_jaccard_cutoff_tp_fp_rates:
    input:
        sam="outputs/cds_fasta_paladin/{acc}_cds.sam",
        csv=expand("outputs/orpheum/{{orpheum_db}}/{{alphabet}}-k{{ksize}}/{{acc}}_{seq}.coding_scores.csv", seq = SEQ)
    output: tsv="outputs/orpheum_cutoffs/{orpheum_db}/{alphabet}-k{ksize}/{acc}.tsv"
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 1
    conda: "envs/tidyverse.yml"
    script: "scripts/compare_jaccard_cutoffs.R"

rule assess_noncoding_in_pseudogenes:
    input: 
        gff = "inputs/gffs/{acc_w_gff}_genomic.gff.gz",
        noncds_bed = "outputs/noncds_bed/{acc_w_gff}_noncds.bed",
        noncds_csv = "outputs/{orpheum_db}/{alphabet}-k{ksize}/{acc_w_gff}_noncds.coding_scores.csv"
    output: 
        pseudo_tsv = "outputs/noncds_pseudogenes/{acc_w_gff}_noncds_predicted_as_coding_in_pseudogenes.tsv",
        pseudo_tab = "outputs/noncds_pseudogenes/{acc_w_gff}_noncds_predicted_as_coding_in_pseudogenes_summary.tsv" 
    resources: 
        mem_mb = 4000,
        tmpdir=TMPDIR
    threads: 1
    conda: "envs/rtracklayer.yml"
    script: "scripts/noncoding_in_pseudogenes.R"

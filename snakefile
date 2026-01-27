import os

# CONFIGURE INPUTS --------------------------------------------------- # 

# Directory with raw reads
READS_DIR = "data/raw_reads/"

# Detect paired FASTQ files
SAMPLES, = glob_wildcards(os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"))
if not SAMPLES:
    raise ValueError(f"No FASTQ files found in {READS_DIR}")

# Path to METABOLIC executable
MB_EXECUTABLE = "/projects/p32449/goop_stirrers/METABOLIC_2025-09-02/METABOLIC/METABOLIC-G.pl"

# Path to GTDB database
GTDB = "/projects/p32449/goop_stirrers/gtdbtk_data/release226"

# -------------------------------------------------------------------- # 

rule all:
    input:
        # GTDB-Tk outputs
        "data/gtdbtk_all",
        # METABOLIC-G outputs
        "data/metabolic_g/METABOLIC_result.xlsx",
        # prokka outputs
        expand("data/prokka/{sample}/{sample}.faa", sample=SAMPLES),
        # EcoFoldDB outputs
        # expand("data/EcoFoldDB_snakemake/annotated/{sample}/{sample}.ecofolddb_annotations.txt", sample=SAMPLES)
        # GToTree outputs
        "data/GToTree/isolate_genome.tre"




rule fastqc_raw:
    input:
        r1 = os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz")
    output:
        "data/fastqc_raw/{sample}_R1_001_fastqc.html",
        "data/fastqc_raw/{sample}_R2_001_fastqc.html",
        "data/fastqc_raw/{sample}_R1_001_fastqc.zip",
        "data/fastqc_raw/{sample}_R2_001_fastqc.zip"
    threads: 4
    shell:
        """
        module load fastqc
        mkdir -p data/fastqc_raw
        fastqc {input.r1} {input.r2} --outdir data/fastqc_raw --threads {threads}
        """


rule multiqc_raw:
    input:
        expand("data/fastqc_raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("data/fastqc_raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES)
    output:
        "data/multiqc_raw/multiqc_report.html"
    shell:
        """
        module load multiqc
        mkdir -p data/multiqc_raw
        multiqc data/fastqc_raw --outdir data/multiqc_raw 
        """

rule trimmomatic:
    input:
        r1 = "data/raw_reads/{sample}_R1_001.fastq.gz",
        r2 = "data/raw_reads/{sample}_R2_001.fastq.gz"
    output:
        r1_paired   = "data/trimmed/{sample}_R1_paired.fastq.gz",
        r1_unpaired = "data/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_paired   = "data/trimmed/{sample}_R2_paired.fastq.gz",
        r2_unpaired = "data/trimmed/{sample}_R2_unpaired.fastq.gz"
    params:
        adapters = "/projects/p31618/databases/adapters.fa"
    threads: 8
    conda:
        "envs/trimmomatic.yml"
    shell:
        """
        mkdir -p data/trimmed
        trimmomatic PE \
          -threads {threads} \
          -Xmx28G \
          {input.r1} {input.r2} \
          {output.r1_paired} {output.r1_unpaired} \
          {output.r2_paired} {output.r2_unpaired} \
          HEADCROP:15 \
          ILLUMINACLIP:{params.adapters}:2:30:10 \
          LEADING:20 TRAILING:20 \
          SLIDINGWINDOW:4:20 MINLEN:50
        """

rule fastqc_trimmed:
    input:
        r1 = "data/trimmed/{sample}_R1_paired.fastq.gz",
        r2 = "data/trimmed/{sample}_R2_paired.fastq.gz"
    output:
        "data/fastqc_trimmed/{sample}_R1_paired_fastqc.html",
        "data/fastqc_trimmed/{sample}_R2_paired_fastqc.html",
        "data/fastqc_trimmed/{sample}_R1_paired_fastqc.zip",
        "data/fastqc_trimmed/{sample}_R2_paired_fastqc.zip"
    threads: 4
    shell:
        """
        module load fastqc
        mkdir -p data/fastqc_trimmed
        fastqc {input.r1} {input.r2} --outdir data/fastqc_trimmed --threads {threads}
        """

rule multiqc_trimmed:
    input:
        expand("data/fastqc_trimmed/{sample}_R1_paired_fastqc.zip", sample=SAMPLES),
        expand("data/fastqc_trimmed/{sample}_R2_paired_fastqc.zip", sample=SAMPLES)
    output:
        "data/multiqc_trimmed/multiqc_report.html"
    shell:
        """
        module load multiqc
        mkdir -p data/multiqc_trimmed
        multiqc data/fastqc_trimmed --outdir data/multiqc_trimmed 
        """

rule spades_wgs:
    input:
        r1 = "data/trimmed/{sample}_R1_paired.fastq.gz",
        r2 = "data/trimmed/{sample}_R2_paired.fastq.gz"
    output:
        contigs = "data/spades/{sample}/contigs.fasta",
        scaffolds = "data/spades/{sample}/scaffolds.fasta"
    threads: 16
    resources:
        slurm_account = "b1042",
        slurm_partition = "genomics-himem",
        runtime = 7200,      
        nodes = 1,
        mem_mb = 150000,       
        slurm_extra = "--mail-user=esmee@u.northwestern.edu --mail-type=END"
    conda:
        "envs/spades.yml"
    shell:
        """
        mkdir -p data/spades/{wildcards.sample}
        spades.py \
          -t {threads} \
          -m 140 \
          -k 21,33,55,77,99 \
          --only-assembler \
          -1 {input.r1} \
          -2 {input.r2} \
          -o data/spades/{wildcards.sample}
        """

rule list_scaffolds:
    input: 
        expand("data/spades/{sample}/scaffolds.fasta", sample=SAMPLES)
    output:
        "data/gtdbtk_input/genome_list.tsv"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as f:
            for sample in SAMPLES:
                genome_path = f"data/spades/{sample}/scaffolds.fasta"
                f.write(f"{genome_path}\t{sample}\n")

rule gtdbtk_classify:
    input:
        "data/gtdbtk_input/genome_list.tsv"
    output:
        directory("data/gtdbtk_all")
    threads: 32
    resources:
        slurm_account="b1042",
        slurm_partition="genomics-himem",
        runtime=48*60,
        nodes=1,
        mem_mb=256000,
        slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
    shell:
        """
        module purge all
        module load python-miniconda3
        source activate /projects/p32449/goop_stirrers/miniconda3/envs/gtdbtk-2.5.2

        gtdbtk classify_wf \
            --batchfile {input} \
            --out_dir {output} \
            --extension scaffolds.fasta \
            --cpus {threads}
        """

rule metabolic_g: 
    input: 
        expand("data/spades/{sample}/scaffolds.fasta", sample=SAMPLES)
    output: 
        directory("data/metabolic_g/METABOLIC_result.xlsx")
    #conda: 
    #    "envs/metabolic_adjusted.yml"
    threads: 8
    resources:
        slurm_account="p32449",
        slurm_partition="short",
        runtime=4*60,
        nodes=1,
        mem_mb=60000,
        slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
    shell: 
        """
        module load perl
        module load python-miniconda3
        
        # Bypass buggy activation script
        export PATH=/projects/p32449/goop_stirrers/miniconda3/envs/METABOLIC_v4.0/bin:$PATH
        export ADDR2LINE=/usr/bin/addr2line

        mkdir -p {output}/scaffolds
        for genome in {input}; do
            folder=$(basename $(dirname "$genome"))
            cp "$genome" "{output}/scaffolds/${{folder}}_scaffolds.fasta"
        done

        export GTDB_DATA_PATH={GTDB}

        perl {MB_EXECUTABLE} \
            -in-gn {output}/scaffolds \
            -o {output} \
            -t {threads}
        """

rule prokka:
    input:
        "data/spades/{sample}/scaffolds.fasta"
    output:
        directory("data/prokka/{sample}/{sample}.faa")
    resources:
        slurm_account="p32449",
        slurm_partition="short",
        runtime=4*60,
        nodes=1,
        mem_mb=60000,
        slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
    shell:
        """
        module load python-miniconda3
        source $(conda info --base)/etc/profile.d/conda.sh
        source activate /projects/p31618/software/prokka
        
        prokka \
            --cpus 2 \
            --outdir {output} \
            --prefix {wildcards.sample} \
            {input}
        """

rule prodigal:
    input:
        "data/spades/{sample}/scaffolds.fasta"
    output:
        a="data/prodigal/{sample}/{sample}.faa",
        o="data/prodigal/{sample}/{sample}.genes",
        d="data/prodigal/{sample}/{sample}.fna"
    resources:
        slurm_account="p32449",
        slurm_partition="normal",
        runtime=24*60, 
        nodes=1,
        mem_mb=16000,
        slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
    shell:
        """
        module load prodigal

        mkdir -p data/prodigal/{wildcards.sample}

        prodigal -i {input} \
             -o {output.o} \
             -a {output.a} \
             -d {output.d} \
             -f gff \
             -p meta
        """

# rule EcoFoldDB:
#     input:
#         "data/prodigal/{sample}/{sample}.faa"
#     output:
#         "data/EcoFoldDB_snakemake/annotated/{sample}/{sample}.ecofolddb_annotations.txt"
#     resources:
#         slurm_account="p32449",
#         slurm_partition="gengpu",
#         runtime=60,
#         nodes=1,
#         mem_mb=10000,
#         gres="gpu:a100:1",
#         slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
#     shell:
#         """
#         module load cuda   

#         # Add Foldseek/EcoFoldDB to PATH
#         export PATH=/projects/p31618/software/EcoFoldDB/foldseek/bin/:$PATH

#         # Run annotation
#         /projects/p31618/software/EcoFoldDB/EcoFoldDB-annotate.sh \
#         --EcoFoldDB_dir /projects/p31618/software/EcoFoldDB/EcoFoldDB_v2.0 \
#         --gpu 1 \
#         --ProstT5_dir /projects/p31618/software/EcoFoldDB/ProstT5_dir \
#         -o data/EcoFoldDB/{wildcards.sample} \
#         {input}
#         """

rule build_tree:
    input:
        ref_genomes="/projects/p32449/maca_mags_metabolic/data/2025-09-16_metagenomes/genome_list_genbank.txt",
        isolate_genomes=expand("data/spades/{sample}/scaffolds.fasta", sample=SAMPLES)
    output:
        "data/GToTree/isolate_genome.tre"
    resources:
        slurm_account="p32449",
        slurm_partition="normal",
        runtime=4*60, 
        nodes=1,
        mem_mb=16000,
        slurm_extra="--mail-user=esmee@u.northwestern.edu --mail-type=END,FAIL"
    shell:
        """
        mkdir -p data/GToTree

        combined_list="data/GToTree/combined_genome_list.txt"
        : > "$combined_list"

        # Convert reference list entries to absolute paths
        while read -r line; do
            [ -z "$line" ] && continue
            realpath "$line" >> "$combined_list"
        done < {input.ref_genomes}

        # Convert isolate genomes to absolute paths
        for g in {input.isolate_genomes}; do
            realpath "$g" >> "$combined_list"
        done

        module load python-miniconda3
        eval "$(conda shell.bash hook)"
        conda activate /projects/p32449/goop_stirrers/miniconda3/envs/gtotree

        GToTree \
            -f "$combined_list" \
            -H Bacteria_and_Archaea \
            -o data/GToTree
        """
include: "Common.snake"

configfile: "config.yaml"


rule all: 
   input:expand("annotation/{assembler}/{group}/{group}_C10K_Cogs_filtered.tsv", group=GROUPS)

#find for ORFs
rule prodigal:
    input:
        "assembly/{assembler}/{group}/{group}.fasta"
    output:
        faa="annotation/{assembler}/{group}/{group}.faa",
        fna="annotation/{assembler}/{group}/{group}.fna",
        gff="annotation/{assembler}/{group}/{group}.gff"
    log: 
        "annotation/{assembler}/{group}/prodigal.log" 
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.fna} -f gff -p meta -o {output.gff} &> {log} "

#cut contigs by ORFs 
rule cut_contigs:
    input:
        fa="assembly/{assembler}/{group}/{group}.fasta",
        gff="annotation/{assembler}/{group}/{group}.gff"
    output:
        contig="annotation/{assembler}/split_c10k/{group}.fasta",
        gff="annotation/{assembler}/{group}/{group}_C10K.gff",
        fna="annotation/{assembler}/{group}/{group}_C10K.fna",
        faa="annotation/{assembler}/{group}/{group}_C10K.faa"
    shell:
        "Use_orf_to_cut.py {input.fa} {input.gff} > {output.contig}"

#Cut faa file in chunks, so that we can have faster annotation  
rule split_fasta:
    input:
        "annotation/{assembler}/{group}/{group}_C10K.faa"
    output:
        touch("annotation/{assembler}/{group}/temp_splits/folder.done")
    threads:
        THREADS
    shell:
        "Fasta_Batchs.py {input} {threads} -t annotation/{wildcards.assembler}/{wildcards.group}/temp_splits/" 


#Do rpsblast annotation  
rule rpsblast_on_folder:
    input:
        name="annotation/{assembler}/{group}/{filename}.faa",
        touch="annotation/{assembler}/{group}/temp_splits/folder.done"
    output:
        "annotation/{assembler}/{group}/{filename}_Rpsblast_cogs.tsv"
    log:
        "annotation/{assembler}/{group}/{filename}_Rpsblast.log"
    threads:
        THREADS
    shell:
        """
        Rpsblast_loop.sh annotation/{wildcards.assembler}/{wildcards.group}/temp_splits {threads} &>log
        mv annotation/{wildcards.assembler}/{wildcards.group}/temp_splits/Complete_Rpsblast_cogs.tsv {output}
        """

# select best hit and use criterion : min 5% coverage, mine 1e-10 evalue
rule parse_cogs_annotation:
    input:
        "annotation/{assembler}/{group}/{filename}_Rpsblast_cogs.tsv"    
    output:
        "annotation/{assembler}/{group}/{filename}_Cogs_filtered.tsv"
    shell:
        "Filter_Cogs.py {input} --cdd_cog_file /home/sebr/seb/Applications/CONCOCT/scgs/cdd_to_cog.tsv  > {output}"









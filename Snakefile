configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001_val_1.fq.gz")

QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',sample=SAMPLES))
print(TRIMMED_FQ)

rule all:
    input:
        QC_FILES,
        TRIMMED_FQ 

rule perform_qc:
    input:
        R1=config['fastq']+'/{sample}_R1_001_val_1.fq.gz',
        R2=config['fastq']+'/{sample}_R2_001_val_2.fq.gz'
    params:
        out_dir = config['output_dir']+'/qc'
    output:
        config['output_dir']+'/qc/{sample}_R1_001_fastqc.html',
        config['output_dir']+'/qc/{sample}_R1_001_fastqc.zip',
        config['output_dir']+'/qc/{sample}_R2_001_fastqc.html',
        config['output_dir']+'/qc/{sample}_R2_001_fastqc.zip'
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input.R1} {input.R2}
        '''

rule trim_galore:
    input:
        R1=config['fastq']+'/{sample}_R1_001_val_1.fq.gz',
        R2=config['fastq']+'/{sample}_R2_001_val_2.fq.gz'
    params:
        out_dir = config['output_dir']+'/trimmed'        
        p = config['nthread']
    output:
        config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',
        config['output_dir']+'/trimmed/{sample}_2_val_2.fq.gz'        
    shell:
        r"""
            trim_galore --paired --cores {params.p} {input} {output} 
        """

# https://stackoverflow.com/questions/46066571/accepting-slightly-different-inputs-to-snakemake-rule-fq-vs-fq-gz

rule align:
    input:
        read = config['data']+"/{sample}.fastq.gz", 
        genome = directory('rnaseq/genome/index/')
    params:
        prefix = 'rnaseq/star/{sample}_'
    output:
        'rnaseq/star/{sample}_Aligned.sortedByCoord.out.bam',
        'rnaseq/star/{sample}_Log.final.out'
    message:
        'mapping {wildcards.sample} to genome'
    shell:
        "mkdir -p {params.outdir}; "
        "cd {params.outdir}; "
        "STAR --runThreadN 4 --genomeDir {input.genome} --readFilesIn {input.read} --readFilesCommand gunzip -c --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --alignEndsType EndToEnd"

# rule make_files:
#     input:
#         dirs = directory('rnaseq/star')
#     output:
#         b1,
#         b2
#     script:
#         "group_files.py"
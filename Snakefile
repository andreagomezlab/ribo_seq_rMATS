configfile:
    "config.json"

SAMPLES_CTRL, = glob_wildcards(config['ctrl']+"/{id}_R1_001_val_1.fq.gz")
print(SAMPLES_CTRL)
SAMPLES_TREAT, = glob_wildcards(config['treatement']+"/{id}_R1_001_val_1.fq.gz")
print(SAMPLES_TREAT)

#EVENTS = ["A3SS", "A5SS", "MXE", "RI", "SE"]. #used 
#JCS = ["JC", "JCEC"] #used for rMats output
#AS_FILES = expand("rnaseq/rmats2/{event}.MATS.{jc}.txt", event = EVENTS, jc = JCS) #rMats output files

#rule all:
#    input:
#        AS_FILES,
#        BAM_FILES,
#        b1,
#        b2

rule all:
    input:        
        expand(config['output_dir']+'/qc/{sample}_R1_001_fastqc.html', sample=SAMPLES_CTRL),
        expand(config['output_dir']+'/qc/{sample}_R2_001_fastqc.html', sample=SAMPLES_CTRL)

rule perform_qc_ctrl:
    input:
        R1 = config['ctrl']+"/{sample}_R1_001_val_1.fq.gz",
        R2 = config['ctrl']+"/{sample}_R2_001_val_2.fq.gz"
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

rule perform_qc_treat:
    input:
        R1 = config['treatement']+"/{sample}_R1_001_val_1.fq.gz",
        R2 = config['treatement']+"/{sample}_R2_001_val_2.fq.gz"
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


# rule align:
#     input:
#         read = config['data']+"/{sample}.fastq.gz", 
#         genome = directory('rnaseq/genome/index/')
#     params:
#         prefix = 'rnaseq/star/{sample}_'
#     output:
#         'rnaseq/star/{sample}_Aligned.sortedByCoord.out.bam',
#         'rnaseq/star/{sample}_Log.final.out'
#     message:
#         'mapping {wildcards.sample} to genome'
#     shell:
#         "mkdir -p {params.outdir}; "
#         "cd {params.outdir}; "
#         "STAR --runThreadN 4 --genomeDir {input.genome} --readFilesIn {input.read} --readFilesCommand gunzip -c --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --alignEndsType EndToEnd"

# rule make_files:
#     input:
#         dirs = directory('rnaseq/star')
#     output:
#         b1,
#         b2
#     script:
#         "group_files.py"

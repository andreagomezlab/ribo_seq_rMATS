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
# https://stackoverflow.com/questions/46066571/accepting-slightly-different-inputs-to-snakemake-rule-fq-vs-fq-gz

def trim_galore_input_determination(wildcards):
     potential_file_path_list = []
     # Cycle through both suffix possibilities:
     for fastqSuffix in [".fq", ".fq.gz"]:

         # Cycle through both read directions
         for readDirection in ['.1','.2']:

             #Build the list for ech suffix
             potential_file_path = config["fq_in_path"] + "/" + wildcards.sample + readDirection + fastqSuffix

             #Check if this file actually exists
             if exists(potential_file_path):

                 #If file is legit, add to list of acceptable files
                 potential_file_path_list.append(potential_file_path)

     # Checking for an empty list
     if len(potential_file_path_list):
         return potential_file_path_list
     else:
         return ["trim_galore_input_determination_FAILURE" + wildcards.sample]

rule trim_galore_unzipped_PE:
    input:
        unpack(trim_galore_input_determination)
    output:
        expand("{trim_out_path}/{{sample}}.{readDirection}.fq.gz",
            trim_out_path=config["trim_out_path"],
            readDirection=['1','2'])
    params:
        out_path=config['trim_out_path'],
    conda:
        'envs/biotools.yaml',
    shell:
        'trim_galore --gzip -o {params.out_path} --paired {input}'



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

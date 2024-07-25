configfile:
    "config.json"

SAMPLES, = glob_wildcards(config['ctrl']+"/{id}_R1_001_val_1.fq.gz")
print(SAMPLES)
EVENTS = ["A3SS", "A5SS", "MXE", "RI", "SE"]. #used 
JCS = ["JC", "JCEC"] #used for rMats output
AS_FILES = expand("rnaseq/rmats2/{event}.MATS.{jc}.txt", event = EVENTS, jc = JCS) #rMats output files

rule all:
    input:
        AS_FILES,
        BAM_FILES,
        b1,
        b2


rule perform_qc:
    input:
        R1=RAWDATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAWDATA_DIR+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir = 'qc'
    output:
       'qc/{sample}_R1_001_fastqc.html',
       'qc/{sample}_R1_001_fastqc.zip',
       'qc/{sample}_R2_001_fastqc.html',
       'qc/{sample}_R2_001_fastqc.zip',
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input.R1} {input.R2}
        '''



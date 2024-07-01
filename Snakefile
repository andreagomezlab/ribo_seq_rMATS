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
# rule aggregate_input_ctrl:
#     input:
#         expand(config['ctrl']+"{sample}_R1_001_val_1.fq.gz:"+config['ctrl']+"{sample}_R2_001_val_2.fq.gz,", sample=SAMPLES, allow_missing=True)       
#     output:
#         config['ctrl']+"/ctrl_input.txt"
#     shell:
#         """
#         echo {input} > {output}
#         """

rule analysesplicing:
    input: 
        r1 = expand('{sample}_R1_001_val_1.fq.gz', sample=SAMPLES),
        r2 = expand('{sample}_R2_001_val_2.fq.gz', sample=SAMPLES)
    params:        
        rlen = config['readLength'], 
        nthread = config['nthread'],
        gtf = config['genome_gtf'],
        star_index = config['star_index'],
        novelSS = config['novel'],
        outdir = config['output_dir']+"/alternative_splicing_rMATS/"+config['project_title'],
        r1 = ','.join(input.r1),
        r2 = ','.join(input.r2)
     conda:
         "binfo-env.yaml"        
    output: config['genome_gtf']+'/alternative_splicing_rMATS/All_T1_vs_T0/MATS_output/AS_Event.MATS.JunctionCountOnly.txt'
    run:        
        shell(r'''source activate python2 && python /home/cmb-panasas2/skchoudh/software_frozen/rMATS.3.2.5/RNASeq-MATS.py -s1 {b1} -s2 {b2} -t paired -len 138 -gtf {gtf} -bi {star_index} -o {params.outdir} 1>&2''')


# python rmats.py --s1 {r1} --s2 {r2} --gtf {gtf} --bi {star_index} -t {paired} --readLength {rlen} --nthread {nthread} --{novelSS} --od $OUT --tmp $OUT_TMP

# rule aggregate_input_treatment:
#     output:
#         "treatment_input.txt"
#     input:
#         expand("{sample}_L004_R1_001_val_1.fq.gz", sample=config['treatement'])
#     shell:
#         """
#         echo {input} > {output}
#         """


# rule make_summary_table:
#     output:
#         "summary_table.txt",
#     input:
#         aggregate_input,
#     shell:
#         """
#         echo {input} >> {output}
#         echo " " >> {output}
#         """

# SAMPLES, = glob_wildcards(config['data']+"/{id}_L004_R1_001_val_1.fq.gz")

# rule move:
#     input:
#         bam = expand("{sample}.final.bam", sample=SAMPLES),
#         bai = expand("{sample}.final.bai", sample=SAMPLES)
#     params:
#         output = config['data']+"/Alignment"
#     shell:
#         """
#         mkdir -p {params.output}
#         mv {input} {params.output}
#         """

# rule align_sort:
#     input:
#         config['data']+"/{sample}_L001_R1_001.fastq.gz",
#         config['data']+"/{sample}_L001_R2_001.fastq.gz"
#     output:
#         temp("{sample}.mapped.bam")
#     params:
#         rg = "@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}",
#         ref = config['genome']
#     conda:
#         "envs.yaml"
#     shell:
#         "bwa mem -R '{params.rg}' -M {params.ref} {input} | samtools sort -o {output} -"

# rule mark_duplicates:
#     input:
#         "{sample}.mapped.bam"
#     output:
#         bam = "{sample}.final.bam",c
#         bai = "{sample}.final.bai",
#         metrics = temp("{sample}_dedup_metricx.txt")
#     conda:
#         "envs.yaml"
#     shell:
#         "picard MarkDuplicates I={input} O={output.bam} M={output.metrics} CREATE_INDEX=true"

# rule analysesplicing:
#     input:
#         b1 = 'rnaseq/star/wt.txt',
#         b2 = 'rnaseq/star/ko.txt',
#         gtf = 'rnaseq/genome/mm.gtf'
#     output:
#         AS_FILES
#     params:
#         outdir = directory('rnaseq/rmats2'),
#         tmp = directory('rnaseq/rmats2/tmp')
#     shell:
#        "mkdir -p {params.outdir}; "
#        "mkdir -p {params.tmp}; "
#        "rmats.py --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --variable-read-length --readLength 100 --libType fr-firststrand --od {params.outdir} --tmp {params.tmp}"
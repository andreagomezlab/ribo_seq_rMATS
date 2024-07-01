configfile:
    "config.json"

print(f'{config["ctrl"]}/')
print(f"{config['ctrl']}/{{id}}fq.gz")

SAMPLES, = glob_wildcards(f"{config['ctrl']}/{{id}}_R1_001_val_1.fq.gz")
print(SAMPLES)


rule aggregate_input_ctrl:
    output:
        "ctrl_input.txt"
    input:
        expand("{sample}_R1_001_val_1.fq.gz", sample=SAMPLES)
    shell:
        """
        echo {input} > {output}
        """

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
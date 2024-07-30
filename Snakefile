configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001_val_1.fq.gz")

QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_2_val_2.fq.gz',sample=SAMPLES))
BAMS=expand(config['output_dir']+'mapped/bams/{sample}.bam', sample=SAMPLES)

print(TRIMMED_FQ)

rule all:
    input:
        QC_FILES,
        TRIMMED_FQ,
        BAMS

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

rule perform_trim_galore:
    input:
        R1=config['fastq']+'/{sample}_R1_001_val_1.fq.gz',
        R2=config['fastq']+'/{sample}_R2_001_val_2.fq.gz'
    params:
        out_dir = config['output_dir']+'/trimmed',        
        p = config['nthread']
    output:
        config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',
        config['output_dir']+'/trimmed/{sample}_2_val_2.fq.gz'        
    shell:
        r"""
            trim_galore --paired --cores {params.p} {input} {output} 
        """

rule perform_STAR_aligner:
    input:
        R1=config['output_dir']+'/trimmed/{sample}_1_val_1.fq.gz',
        R2=config['output_dir']+'/trimmed/{sample}_2_val_2.fq.gz',     
        index=config['star_index']
    output: config['output_dir']+'mapped/bams/{sample}.bam'
    params:
        prefix = config['output_dir']+'mapped/bams/{sample}',
        unmapped = config['output_dir']+'unmapped/fastq/{sample}',
        starlogs = config['output_dir']+'mapped/starlogs'
    threads: config['nthread']
    shell:
        r'''
        STAR --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --outFilterMatchNmin 50\
             --outFilterMismatchNmax 100\
             --readFilesCommand zcat\
             --outReadsUnmapped {params.unmapped} && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
        '''

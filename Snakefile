configfile:
    "config.json"

(SAMPLES,) = glob_wildcards(config['fastq']+"/{id}_R1_001.fastq.gz")

EVENTS = ["A3SS", "A5SS", "MXE", "RI", "SE"] #used for rMats output
JCS = ["JC", "JCEC"] #used for rMats output
b1 = config['output_dir']+"/rMATs/"+config['project_title']+"/ctrl.txt"
b2 = config['output_dir']+"/rMATs/"+config['project_title']+"/treat.txt"
QC_FILES=expand(config['output_dir']+"/qc/{sample}_R1_001_fastqc.html", sample=SAMPLES)
TRIMMED_FQ=expand(config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',sample=SAMPLES)
TRIMMED_FQ.append(expand(config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',sample=SAMPLES))
BAMS=expand(config['output_dir']+'/mapped/bams/{sample}.bam', sample=SAMPLES)
print(BAMS)
AS_FILES = expand(config['output_dir']+"/rMATs/"+config['project_title']+"/{event}.MATS.{jc}.txt", event = EVENTS, jc = JCS) #rMats output files

rule all:
    input:
        QC_FILES,
        TRIMMED_FQ,
        BAMS,
       	b1,
        b2,
        AS_FILES

rule merge_files_R1:
    input:
        F1=config['fastq']+'/{sample}_R1_001.fastq.gz',
        F2=config['fastq']+'/{sample}_1.fq.gz'
    params:
        out_dir = config['output_dir']+'/merged'
    output:
        config['output_dir']+'/merged/{sample}_merged_2.fastq.gz'
    shell :
        r'''
            cat {input.F1} {input.F2} > {output}
        '''    

rule merge_files_R2:
    input:
        F1=config['fastq']+'/{sample}_R2_001.fastq.gz',
        F2=config['fastq']+'/{sample}_2.fq.gz'
    params:
        out_dir = config['output_dir']+'/merged'
    output:
        config['output_dir']+'/merged/{sample}_merged_2.fastq.gz'
    shell :
        r'''
            cat {input.F1} {input.F2} > {output}
        '''    

rule perform_qc:
    input:
        R1=config['fastq']+'/{sample}_R1_001.fastq.gz',
        R2=config['fastq']+'/{sample}_R2_001.fastq.gz'
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
        R1=config['fastq']+'/{sample}_R1_001.fastq.gz',
        R2=config['fastq']+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir = config['output_dir']+'/trimmed'                
    output:
        config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',                    
        config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz'        
    shell:
        r"""
            trim_galore --paired --phred33 --cores 8 {input.R1} {input.R2} -o {params.out_dir} 
        """

rule perform_STAR_aligner:
    input:
        R1=config['output_dir']+'/trimmed/{sample}_R1_001_val_1.fq.gz',
        R2=config['output_dir']+'/trimmed/{sample}_R2_001_val_2.fq.gz',     
        index=config['star_index']
    output: 
        config['output_dir']+'/mapped/bams/{sample}.bam'
    params:
        prefix = config['output_dir']+'/mapped/bams/{sample}',
        unmapped = config['output_dir']+'/unmapped/fastq/{sample}',
        starlogs = config['output_dir']+'/mapped/starlogs'
    threads: 
        config['nthread']
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

rule make_files:
    input:                
        bam_path = config['output_dir'] +'/mapped/bams/',
        metadata = config['metadata']
    output:
        b1,
        b2
    script:
        "scripts/group_input_rMATs.py"


rule run_rMATS:
    input:
        b1 = config['output_dir']+"/rMATs/"+config['project_title']+"/ctrl.txt",
        b2 = config['output_dir']+"/rMATs/"+config['project_title']+"/treat.txt",
        gtf = config['genome_gtf']
    output:
        AS_FILES
    params:
        outdir = config['output_dir']+"/rMATs/"+config['project_title'],
        tmp = config['output_dir']+"/rMATs/"+config['project_title']+"/tmp",
        readTy = config['readType'],
        readLen = config['readLength'],
        novel = config['novel'],
        nt = config['nthread']
    conda:
        config['rMATs_environment']
    shell:       
       "mkdir -p {params.tmp}; "
       "python /apps/rmats-turbo/rmats.py --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t {params.readTy} --readLength {params.readLen} --nthread {params.nt} --{params.novel} --od {params.outdir} --tmp {params.tmp}"




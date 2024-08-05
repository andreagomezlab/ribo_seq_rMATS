

## Step 1: Install rmats-turbo

Install [rMATs-turbo](https://github.com/Xinglab/rmats-turbo) using conda environment option and inform the enviroment name in <code>config.json</code>  considering the `rMATs_environment` propertie.

## Step 2: Install ribo_seq_rMATS workflow

### 1. Clone repository
```shell
$ git clone https://github.com/andreagomezlab/ribo_seq_rMATS.git
$ cd ribo_seq_rMATS/
```

### 2. Create and activate the conda environment
```shell
$ cd env/
$ conda env create --name envname --file=scrna_env.yml 
$ conda activate rnaseq_env
```

### 3. Dry run the pipeline workflow
```shell
$ snakemake -np
```

## Step 3: Configure workflow

Configure the workflow if necessary by editing the file <code>config.json</code> and <code>metadata.cvs</code>


## Step 4: Execute the workflow

```shell
$ snakemake -c10 --use-conda
$ RScripts

```

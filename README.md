
## Step 1: Install workflow

### 1. Clone repository
```shell
$ git clone https://github.com/andreagomezlab/ribo_seq_rMATS.git
$ cd ribo_seq_rMATS
```

### 2. Create and activate the conda environment
```shell
$ cd env/
$ conda env create --name envname --file=scrna_env.yml 
$ conda activate rnaseq_env
```

### 3. Dry run the pipeline workflow
```shell
$ snakemake --np
```

## Step 2: Configure workflow

Configure the workflow if necessary by editing the file <code>config.json</code>


## Step 3: Execute the workflow
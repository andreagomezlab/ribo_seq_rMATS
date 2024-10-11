import os
import pandas as pd

def make_files(bam_path, metadata, out1, out2):    
    meta = pd.read_csv(metadata, sep=',', header=None)
    print(meta)
    dictionary = {"ctrl": [], "treat": []} 

    for index, row in meta.iterrows():
        print(row)
        if row[1] == "ctrl":
            dictionary["ctrl"].append(bam_path + '/' + row[0]+".bam")
        else:
            dictionary["treat"].append(bam_path + '/' + row[0]+".bam")

    print(dictionary)

    for key, value in dictionary.items():
        with open (out2 if key == "ctrl" else out1, 'w') as f:
                f.write(','.join(value)) 
                f.close

make_files(snakemake.params[0], snakemake.params[1], snakemake.output[0], snakemake.output[1])


import os

def make_files(data_path, out1, out2):  
    files = []
    substring = "IRS"
    dictionary = {"ko": [], "wt": []} 

    for file in os.listdir(data_path):
        if file.endswith("Aligned.sortedByCoord.out.bam"):
            files.append(file)
     
    for file in files:
        if substring in file:  
            dictionary["ko"].append(file)  
        else:
            dictionary["wt"].append(file) 

    for key, value in dictionary.items():
        with open (out2 if key == "wt" else out1, 'w') as f:
                f.write(','.join(value)) 
                f.close

make_files(snakemake.input[0], snakemake.output[0], snakemake.output[1])

#in the rule that executes the script, I will define the folder with bam files as the input, 
# and final output files as two outputs
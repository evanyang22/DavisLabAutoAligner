#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor
#fasta_file="p22323-s013_EV-Justin-Pre-Op-miRNA_S25_L002_R1_001_trimmed.fq"
#species="oar"
#file_path="p22323-s013_EV-Justin-Pre-Op-miRNA_S25_L002_R1_001_trimmed.fq"
#genome= # unsure
#mapper_output="test.fa"
#mapped_file= #unsure
#precursors= "hairpin.fa"
#matures= "mature.fa"

def filter_fasta(fasta_file, species):
    f=open(fasta_file,'r')
    lines=f.readlines()
    results=[]
    i=0
    while i<len(lines):
        if species in lines[i]:
            results.append(lines[i])
            results.append(lines[i+1])
        i+=1  
    output_file=species + "_"+ os.path.basename(fasta_file)
    ofile = open(output_file, "w")
    for item in results:
        ofile.write(item)
    ofile.close()
    return ofile
def run_mapper(file_path, genome, mapper_output):
    command = ["mapper.pl", file_path, "-e", "-h", "-j", "-m", "-s", mapper_output, "-p", genome, "-t", mapper_output + ".arf"]
    subprocess.run(command, check=True)
def run_quantifier(mapped_file, precursors, matures):
    command = ["quantifier.pl", "-p", precursors, "-m", matures, "-r", mapped_file]
    subprocess.run(command, check=True)
def main():
    parser = argparse.ArgumentParser(description="Run miRDeep2 on a directory of trimmed FastQ files.")
    parser.add_argument("input_dir", help="Directory containing trimmed FastQ files.")
    parser.add_argument("output_dir", help="Directory to write output files.")
    parser.add_argument("genome", help="Genome file for miRDeep2.")
    parser.add_argument("precursors", help="All species precursor miRNA sequences for miRDeep2.")
    parser.add_argument("matures", help="All species mature miRNA sequences for miRDeep2.")
    parser.add_argument("species", help="Species abbreviation for miRDeep2.")
    args = parser.parse_args()
    precursors_species = filter_fasta(args.precursors, args.species)
    matures_species = filter_fasta(args.matures, args.species)
    print(precursors_species)
    print(matures_species)
    files = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir) 
            if file.endswith('_R1_001_trimmed.fq') or file.endswith('_R1_trimmed.fq')]
    print(files)
    mapper_files = [os.path.join(args.output_dir, os.path.basename(file)+"_mapped") for file in files]
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(run_mapper, files, [args.genome]*len(files), mapper_files)
    mapped_files = [file for file in mapper_files]
    print("Created Mapped Files")
    print(mapped_files)
    #with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    #    executor.map(run_quantifier, mapped_files, [precursors_species]*len(mapped_files), [matures_species]*len(mapped_files))
    os.system("quantifier.pl -p " + precursors_species.name + " -m "+ matures_species.name + " -r " +mapped_files[0])
    print("Quantifier Done")
if __name__ == "__main__":
    main()


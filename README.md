# DavisLabAutoAligner
Uses Star to automatically align fastq files

Documentation

**How to Use STAR General Instructions**
1. Install trimgalore, fastqc, star,cutadapt
2. Download sheep reference genome from Ensembl database: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000298735.2/ 
Note: This can be any reference genome, also download the annotation files in GTF format
3. Activate cutadapt env : “conda activate cutadaptenv”
4. Run a quality control : “fastqc input.fastq”
5. Trim adapters and low quality bases from the input fastq file : “trim_galore  --quality 20 --fastqc --illumina input.fastq”
6. Generate reference genome index file (Use GCF ones)

Note : Files with “.fna” may need to be renamed to “.fa”
Note : Better to use full directories rather than abbreviated ones
Command : “STAR --runThreadN 48 -- runMode genomeGenerate --genomeDir genomeDirectory --genomeFastaFiles _directory to fasta file(s)_ --sjdbGTFfile _directory to genome GTF file_" 
Note : Only needs to be done once
7. Use STAR to align genome to reference genome
Command : “STAR --genomeDir _directory to reference genome_ --readFilesIn _directory to fastq files_ --quantMode GeneCounts”

In this case, directory to reference genome is : “/home/ehyang4/ncbi_dataset/data”

**How to Use My Program**
1. Unzip all files into .fastq format
2. Open terminal and cd to afolder with all the fastq files
3. Activate cutadapt env
4. Start Python in terminal via "python" 
5. Copy and paste the _#finished loop program_ from the jupyter notebook into terminal
6. Run it and it should begin processing all the fastq files in the folder and store them into an output folder automatically


Reference Genomes
1. https://www.nature.com/articles/s41437-018-0090-1- says to use OAR 3.1 for miRNA



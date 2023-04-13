# DavisLabAutoAligner
Uses Star to automatically align fastq files

Documentation

**How to Use STAR General Instructions**
1. Install trimgalore, fastqc, star,cutadapt
2. Download sheep reference genome from Ensembl database: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000298735.2/ 
Note: This can be any reference genome
3. Activate cutadapt env : “conda activate cutadaptenv”
4. Run a quality control : “fastqc input.fastq”
5. Trim adapters and low quality bases from the input fastq file : “trim_galore  --quality 20 --fastqc --illumina input.fastq”
6. Generate reference genome index file
Note : Files with “.fna” may need to be renamed to “.fa”
Note : Better to use full directories rather than abbreviated ones
Command : “STAR --runThreadN 48 -- runMode genomeGenerate --genomeDir genomeDirectory --genomeFastaFiles _directory to fasta file(s)_ --sjdbGTFfile _directory to genome GTF file_" 
Note : Only needs to be done once
7. Use STAR to align genome to reference genome
Command : “STAR --genomeDir _directory to reference genome_ --readFilesIn _directory to fastq files_ --quantMode GeneCounts”

In this case, directory to reference genome is : “/home/ehyang4/ncbi_dataset/data”






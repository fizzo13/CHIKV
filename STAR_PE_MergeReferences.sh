#!/bin/bash
#SBATCH --job-name=STAR_PE
#SBATCH --partition=bigmem
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fizzo@nygenome.org
#SBATCH --mem=100g
#SBATCH --time=24:00:00
#SBATCH --output=stdout_%j.log
#SBATCH --cpus-per-task=12


module load star/2.5.2a

# Generate indexed genome
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /gpfs/commons/home/fizzo/MGN/Genomes/mm10_CHIKV_merge --genomeFastaFiles /gpfs/commons/home/fizzo/MGN/Genomes/mm10/mm10.fa /gpfs/commons/home/fizzo/MGN/Genomes/CHIKV/CHIKV.genome.fa --sjdbGTFfile /gpfs/commons/home/fizzo/MGN/Genomes/mm10/gencode.vM25.annotation.gtf

# Read in arguments
STAR_DIR=$1

for files in $(cat /gpfs/commons/home/fizzo/MGN/Scripts/File_index_trimmed.txt); do 
STAR --genomeDir /gpfs/commons/home/fizzo/MGN/Genomes/mm10_CHIKV_merge/ \
--readFilesIn /gpfs/commons/home/fizzo/MGN/TrimmedFASTQ/${files}_L002_R1_001.fastq.gz.trimmed.fa  /gpfs/commons/home/fizzo/MGN/TrimmedFASTQ/${files}_L002_R2_001.fastq.gz.trimmed.fa \
--outFileNamePrefix /gpfs/commons/home/fizzo/MGN/BAM_mergeReference/${files} \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 12 \
--twopassMode Basic; done

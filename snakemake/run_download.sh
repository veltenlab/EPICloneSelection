#$ -q short-sl7
#$ -N download_SRR1037476
#$ -e /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/SRR1037476/fastq/stderr.txt 
#$ -o /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/SRR1037476/fastq/stdout.txt 
#$ -l virtual_free=40G 
#$ -l h_rt=05:55:00 

sample=SRR1037476
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/snakemake
/users/lvelten/mscherer/conda/envs/sra/bin/prefetch ${sample} --output-directory /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/fastq/ --max-size 50G
/users/lvelten/mscherer/conda/envs/sra/bin/fastq-dump --split-3 /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/fastq/${sample}/${sample}.sra -O /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/fastq/


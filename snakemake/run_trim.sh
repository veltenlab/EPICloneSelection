#$ -q short-sl7
#$ -N trim_SRR1037476
#$ -e /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/SRR1037476/trimmed/stderr.txt 
#$ -o /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/SRR1037476/trimmed/stdout.txt 
#$ -l virtual_free=40G 
#$ -l h_rt=05:55:00 

sample=SRR1037476
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/genotools
trim_galore --nextera --paired /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/fastq/${sample}_1.fastq /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/fastq/${sample}_2.fastq -o /nfs/no_backup/lvelten/mscherer/projects/Methylome/data/external/Cabezas/raw/${sample}/trimmed

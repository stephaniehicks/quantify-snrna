#$ -pe local 4
#$ -cwd
#$ -l mem_free=15G,h_vmem=20G

# create salmon index (this process takes ~2-3 hours)
# make sure you are using the right index (for mRNA, pre-mRNA, or both)
# salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/GRCh38.premRNA.fa.gz -i /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0-premRNA --gencode --threads 4
salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/gencode.v32.transcripts.fa.gz -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0-mRNA --gencode --threads 4

#$ -pe local 4
#$ -cwd
#$ -l mem_free=40G,h_vmem=45G

# create salmon index (this process takes ~3 hours)
# salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gencode.v32.preandmrna.fa.gz -i /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0 --threads 4
salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/gencode.v32.preandmrna.fa.gz -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0 --threads 4

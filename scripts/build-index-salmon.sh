#$ -pe local 4
#$ -cwd
#$ -l mem_free=40G,h_vmem=45G

cd /users/shicks1/data/reference/genomes/hsapiens/GENCODE/release_32

# create salmon index (3-5 mins)
salmon index --gencode -t gencode.v32.preandmrna.fa.gz -i /fastscratch/myscratch/shicks1/alsf-filbin/gencode.v32_salmon-index-v1.0.0 --keepDuplicates --threads 4

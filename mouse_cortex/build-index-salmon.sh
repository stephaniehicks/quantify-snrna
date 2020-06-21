#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=10G,h_fsize=30G

# create salmon index (this process takes ~2-3 hours)
salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/gentrome_transcripts.fa.gz \
                -d /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/decoys.txt \
                -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/salmon_index_gentrome_decoys \
                --gencode --threads 4


# default k-mer size is 31; trying a smaller k-mer size
salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/gentrome_transcripts.fa.gz \
                -d /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/decoys.txt \
                -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/salmon_index_gentrome_decoys_k25 \
                --gencode --threads 4 -k 25

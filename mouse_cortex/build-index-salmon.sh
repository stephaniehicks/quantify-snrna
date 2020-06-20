#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=10G,h_fsize=30G

# create salmon index (this process takes ~2-3 hours)
# make sure you are using the right transcriptome (for mRNA, pre-mRNA, intron, or some combination)
# salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/GRCm38.premRNA.fa.gz \
#                 -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/salmon_transcripts_index_mouse_premRNA --gencode --threads 4
salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/gentrome_transcripts.fa.gz \
                -d /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/decoys.txt \
                -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/salmon_transcripts_index --gencode --threads 4

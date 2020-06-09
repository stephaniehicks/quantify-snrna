#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G

# create salmon index (this process takes ~2-3 hours)
# make sure you are using the right transcriptome (for mRNA, pre-mRNA, intron, or some combination)
# salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/GRCh38.premRNA.fa.gz -i /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0-premRNA --gencode --threads 4
salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gentrome_transcripts_mouse.fa.gz \
                -d /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/decoys_mouse.txt \
                -i /fastscratch/myscratch/akuo/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0-transcripts-mouse --gencode --threads 4

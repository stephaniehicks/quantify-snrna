#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=30G,h_vmem=30G,h_fsize=30G

# mem_free=30G,h_vmem=30G,h_fsize=30G (preandmrna index)
# mem_free=25G,h_vmem=25G,h_fsize=30G (intron index)
# create salmon index (this process takes ~2-3 hours for humans)

pipeline=$1 # first parameter (options are: transcripts, preandmran, introncollapse, or intronseparate)
salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/gentrome_${pipeline}_mouse.fa.gz \
             -d /fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/decoys_mouse.txt \
             -i /fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/gencode.vM25_salmon-index-v1.0.0-${pipeline}-mouse \
             --gencode --threads 4


# create salmon index without decoys
# salmon index -t /fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/gencode.vM25.transcripts.fa.gz \
#              -i /fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/salmon_files/gencode.vM25_salmon-index-v1.0.0-transcripts-mouse-nodecoys \
#              --gencode --threads 4

# default k-mer size is 31; trying a smaller k-mer size
# salmon index -t /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/gentrome_transcripts.fa.gz \
#                 -d /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/decoys.txt \
#                 -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/mouse/salmon_index_gentrome_decoys_k25 \
#                 --gencode --threads 4 -k 25

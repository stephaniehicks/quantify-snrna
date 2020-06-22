#$ -pe local 10
#$ -R y
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=10G

# set project directory path
#d=/fastscratch/myscratch/shicks1/alsf-filbin
#f=/fastscratch/myscratch/shicks1/sra
d=/fastscratch/myscratch/akuo/alsf-filbin

# Albert's
salmon alevin --libType ISR \
      --index $d/salmon_files/gencode.v32_salmon-index-v1.0.0-transcripts-mouse \
      -1 $d/SRR9169236.sra_1.fastq \
      -2 $d/SRR9169236.sra_2.fastq \
      --tgMap $d/salmon_files/gencode.v32.annotation.tx2gene.mouse.txt \
      --chromium \
      --threads 10 \
      --output $d/ \
      --dumpFeatures --dumpBfh


# Stephanie's
# salmon alevin --libType ISR \
#       --index $d/salmon_files/mouse/salmon_index_gentrome_decoys_k25 \
#       -1 $f/SRR9169228_1.fastq \
#       -2 $f/SRR9169228_2.fastq \
#       --tgMap $d/salmon_files/mouse/gencode.v32.annotation.tx2gene.mouse.txt \
#       --chromium \
#       --threads 10 \
#       --output $d/mouse_cortex/salmon_quants/ \
#       --dumpFeatures --dumpBfh

# for fn in `cat $d/sample_data/unique_cell_paths.txt`; 
# do 
# samp=`basename ${fn}`
# samp=${samp::-2}
# echo "Processing sample ${samp}"
# salmon quant -i $d/salmon_files/gencode.v32_salmon-index-v1.0.0 -l A \
#          -1 ${fn}1.fastq.gz \
#          -2 ${fn}2.fastq.gz \
#          -p 6 --validateMappings -o $d/salmon_quants/${samp}_quant
# done

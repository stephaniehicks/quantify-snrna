#$ -pe local 2
#$ -R y
#$ -cwd
#$ -l mem_free=60G,h_vmem=60G

# set project directory path
# d=/fastscratch/myscratch/shicks1/alsf-filbin
d=/fastscratch/myscratch/akuo/alsf-filbin

salmon alevin -l ISR -1 $d/SRR9169236.sra_1.fastq -2 $d/SRR9169236.sra_2.fastq --chromium -i $d/salmon_files/gencode.v32_salmon-index-v1.0.0-transcripts-mouse -p 2 -o $d/alevin_output --tgMap $d/salmon_files/gencode.v32.annotation.tx2gene.mouse.txt


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

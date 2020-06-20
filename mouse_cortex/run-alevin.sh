#$ -pe local 10
#$ -R y
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=10G

# set project directory path
d=/fastscratch/myscratch/shicks1/alsf-filbin
f=/fastscratch/myscratch/shicks1/sra
# d=/fastscratch/myscratch/akuo/alsf-filbin

salmon alevin -l ISR \ 
      -1 $f/SRR9169228_1.fastq \ 
      -2 $f/SRR9169228_2.fastq \ 
      --chromium  \ 
      -i $d/salmon_files/mouse/salmon_transcripts_index \ 
      -p 10 \ 
      -o $d/mouse_cortex/salmon_quants/ \ 
      --tgMap $d/salmon_files/mouse/gencode.v32.annotation.tx2gene.mouse.txt 
      


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

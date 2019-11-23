#$ -pe local 6
#$ -cwd
#$ -l mem_free=60G,h_vmem=65G

for fn in `cat /users/shicks1/data/alsf_filbin/sample/unique_cell_filenames.txt | head -2`; 
do 
samp=`basename ${fn}`
samp=${samp::-2}
echo "Processing sample ${samp}"
salmon quant -i /users/shicks1/data/reference/genomes/hsapiens/GENCODE/release_32/gencode.v32_salmon-index-v0.10.2 -l A \
         -1 ${fn}1.fastq.gz \
         -2 ${fn}2.fastq.gz \
         -p 6 --validateMappings -o /fastscratch/myscratch/shicks1/alsf-filbin/salmon_quants/${samp}_quant
done
#$ -pe local 6
#$ -cwd
#$ -l mem_free=60G,h_vmem=65G

for fn in `cat /fastscratch/myscratch/shicks1/alsf-filbin/unique_cell_paths.txt | head -2`; 
do 
samp=`basename ${fn}`
samp=${samp::-2}
echo "Processing sample ${samp}"
salmon quant -i /users/shicks1/data/reference/genomes/hsapiens/ENSEMBL/GRCh38/salmon/ensemble.grch38_salmon-index-v0.14.1 -l A \
         -1 ${fn}1.fastq.gz \
         -2 ${fn}2.fastq.gz \
         -p 6 --validateMappings -o /fastscratch/myscratch/shicks1/alsf-filbin/salmon_quants/${samp}_quant
done
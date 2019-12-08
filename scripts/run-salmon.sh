#$ -pe local 6
#$ -cwd
#$ -l mem_free=40G,h_vmem=45G

for fn in `cat /fastscratch/myscratch/shicks1/alsf-filbin/sample_data/unique_cell_filenames.txt`; 
do 
samp=`basename ${fn}`
samp=${samp::-2}
echo "Processing sample ${samp}"
salmon quant -i /fastscratch/myscratch/shicks1/alsf-filbin/salmon_files/gencode.v32_salmon-index-v1.0.0 -l A \
         -1 ${fn}1.fastq.gz \
         -2 ${fn}2.fastq.gz \
         -p 6 --validateMappings -o /fastscratch/myscratch/shicks1/alsf-filbin/salmon_quants/${samp}_quant
done
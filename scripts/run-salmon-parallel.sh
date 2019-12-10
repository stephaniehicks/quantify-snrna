#$ -pe local 6
#$ -cwd
#$ -l mem_free=10G,h_vmem=15G
#$ -t 1-576
#$ -tc 30

# set project directory path
# d=/fastscratch/myscratch/shicks1/alsf-filbin
d=/fastscratch/myscratch/akuo/alsf-filbin

# run salmon
samplefile=$d/sample_data/unique_cell_paths.txt;
fn=`awk -F'\r' -v var=$SGE_TASK_ID '{if(NR==var)print $1}' $samplefile`;
samp=`basename ${fn}`
samp=${samp::-2}
echo "Processing sample ${samp}"
salmon quant -i $d/salmon_files/gencode.v32_salmon-index-v1.0.0 -l A \
         -1 ${fn}1.fastq.gz \
         -2 ${fn}2.fastq.gz \
         -p 6 --validateMappings -o $d/salmon_quants/${samp}_quant
done
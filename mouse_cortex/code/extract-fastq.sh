#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=5G,h_vmem=10G

a=/fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex
b=/fastscratch/myscratch/shicks1/alsf-filbin


for file in $(<$a/files/SRR_files_10x.txt)
do
    echo "$file.sra"
    # stephanie
    # fastq-dump --split-files --gzip $b/sra/$file
    
    # albert
    fasterq-dump -O $c/sample_data/geo/sra/ \
       -f --include-technical \
       --split-files --gzip $c/sample_data/geo/sra/$file
done
#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=15G,h_fsize=100G

a=/fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex
b=/fastscratch/myscratch/shicks1/alsf-filbin


for file in $(<$a/files/SRR_files_10x.txt)
do
    echo "$file.sra"
    # stephanie
    # fastq-dump --split-files --gzip $b/sra/$file
    
    # albert
    fasterq-dump -O $a/sample_data/geo/sra/ \
       -f --include-technical \
       --split-files $a/sample_data/geo/sra/$file.sra
done

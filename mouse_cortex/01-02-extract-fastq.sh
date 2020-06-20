#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=40G,h_vmem=40G

d=/fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/00-files/

for file in $(<$d/SRR_files_10x.txt)
do
    echo "$file.sra"
     fastq-dump -O "/fastscratch/myscratch/akuo/alsf-filbin/sample_data/geo/sra/" --split-files "/fastscratch/myscratch/akuo/alsf-filbin/sample_data/geo/sra/$file.sra"
done

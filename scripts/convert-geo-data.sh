#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=40G,h_vmem=40G

d=/fastscratch/myscratch/akuo/alsf-filbin/files

for file in $(<$d/SRR_files.txt)
do
    echo "$file.sra"
    fasterq-dump -O "/fastscratch/myscratch/akuo/alsf-filbin/sample_data/geo/sra/" --split-files "/fastscratch/myscratch/akuo/alsf-filbin/sample_data/geo/sra/$file.sra"
done

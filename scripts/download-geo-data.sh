#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G

prefetch --option-file SRR_files.txt
cat SRR_files.txt | parallel -j 4 fasterq-dump {}.sra

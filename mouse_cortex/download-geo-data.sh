#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G

d=/fastscratch/myscratch/akuo/alsf-filbin/mouse_cortex/00-files/

prefetch --option-file $d/SRR_files_10x.txt

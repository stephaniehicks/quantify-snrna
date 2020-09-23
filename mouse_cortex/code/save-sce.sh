#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=40G,h_vmem=40G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R
R CMD BATCH save-sce.R

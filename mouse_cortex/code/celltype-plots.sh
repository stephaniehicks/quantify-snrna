#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=20G,h_vmem=20G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R/4.0
R CMD BATCH celltype-plots.R

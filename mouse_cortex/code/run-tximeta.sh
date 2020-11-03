#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=30G,h_vmem=30G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R/4.0
R CMD BATCH run-tximeta.R

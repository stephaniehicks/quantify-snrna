#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=15G,h_vmem=15G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R/4.0
R CMD BATCH cqn.R

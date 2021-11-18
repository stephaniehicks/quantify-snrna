#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G
#$ -t 1:56

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R/4.0
R CMD BATCH distribution-plots.R

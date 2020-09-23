#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=35G,h_vmem=35G
#$ -t 1:56

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R
R CMD BATCH distribution-plots.R

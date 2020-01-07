#$ -R y
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G
#$ -t 1-534
#$ -tc 50

# set project directory path
# d=/fastscratch/myscratch/shicks1/alsf-filbin
d=/fastscratch/myscratch/akuo/alsf-filbin

# make log directory for error and output files
mkdir -p $d/scripts/log

# run script
module load conda_R
R CMD BATCH compute-mle.R

#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=2G,h_vmem=2G
#$ -t 1:31


# run script
module load conda_R/4.0
R CMD BATCH scrna-distribution-sims.R

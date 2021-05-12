#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=8G,h_vmem=8G
#$ -t 1:31


# run script
module load conda_R/4.0
R CMD BATCH scrna-distribution-sims.R

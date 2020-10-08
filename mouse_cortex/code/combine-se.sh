#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R
R CMD BATCH combine-se.R

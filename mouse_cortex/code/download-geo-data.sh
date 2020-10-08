#$ -cwd
#$ -o log/
#$ -e log/

module load conda_R
R CMD BATCH download-geo-data.R

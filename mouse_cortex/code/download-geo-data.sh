#$ -cwd
#$ -o log/
#$ -e log/

module load conda_R/4.0
R CMD BATCH download-geo-data.R

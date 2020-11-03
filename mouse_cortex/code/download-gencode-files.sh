#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=40G,h_vmem=40G

# set project directory path
d=/fastscratch/myscratch/akuo/alsf-filbin

# run script
module load conda_R/4.0
R CMD BATCH download-gencode-files.R

# create decoy-aware transcriptome
cd ../salmon_files/
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys_mouse.txt
sed -i.bak -e 's/>//g' decoys_mouse.txt

cat gencode.vM25.transcripts.mouse.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome_transcripts_mouse.fa.gz
cat gencode.vM25.preandmrna.mouse.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome_preandmrna_mouse.fa.gz
cat gencode.vM25.introncollapse.mouse.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome_introncollapse_mouse.fa.gz
cat gencode.vM25.intronseparate.mouse.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome_intronseparate_mouse.fa.gz

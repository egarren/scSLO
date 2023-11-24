##### Download data ######

#### Cell ranger pipeline (for 10X experiments)
module load cellranger

cd //GSE203132

sed -i -e 's/\r$//' GSE203132.txt
for d in $(cat GSE203132.txt); do
     cellranger count --id=$d --nosecondary --transcriptome=/ref_genomes/mm39 \
--fastqs=/RNAseq_fastq/GSE203132/$d
done


sed -i -e 's/\r$//' GSE203132_2.txt
for d in $(cat GSE203132_2.txt); do
     cellranger count --id=$d --nosecondary --transcriptome=/ref_genomes/mm39 \
--fastqs=/RNAseq_fastq/GSE203132/$d,/RNAseq_fastq/GSE203132/"$d"_2
done

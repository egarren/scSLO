##install velocyto
pip3 install velocyto
cd /ref_genomes/mm39/genes/
gunzip genes.gtf.gz

cd /RNAseq_cellranger/B/GSE203132
sed -i -e 's/\r$//' GSE203132_3.txt
for d in $(cat GSE203132_3.txt); do
     velocyto run10x --samtools-threads 8 --samtools-memory 28000 -m /ref_genomes/mm39_rmsk.gtf /scB/cellranger2/$d /ref_genomes/mm39/genes/genes.gtf
done

cd /scB/
find cellranger2/ -type f | grep -i loom$ | xargs -i cp {} loom


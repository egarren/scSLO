cp -r /RNAseq_cellranger/scTfh/gex_2020 /scTFH/cellranger

cd /RNAseq_cellranger/scTfh/gex_2020 
sed -i -e 's/\r$//' samples.txt
for d in $(cat samples.txt); do
      velocyto run10x --samtools-threads 8 --samtools-memory 28000 -m /ref_genomes/mm39_rmsk.gtf /scTFH/cellranger/$d /ref_genomes/mm39/genes/genes.gtf
done

cd /scTFH
find cellranger/ -type f | grep -i loom$ | xargs -i cp {} loom

rm -r cellranger

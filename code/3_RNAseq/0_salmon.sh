
#create genome index
mkdir GRCm39
cd ./GRCm39
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip GCF_000001635.27_GRCm39_genomic.gtf.gz
mkdir ref
cd ./ref
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna --sjdbGTFfile GCF_000001635.27_GRCm39_genomic.gtf --sjdbOverhang 99

#salmon index
cd /FDC_rnaseq
mkdir salmon
cd salmon
wget http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/salmon_sa_index?tag=default
tar zxvf salmon_sa_index.tgz
mv default mm10_index
salmon index -t ../GRCm39/GCF_000001635.27_GRCm39_genomic.fna -i mm39_index

#salmon count
cd /FDC_rnaseq
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
    do salmon quant -i ../../../../salmon/mm10_index -l A -r "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz "$i"_L003_R1_001.fastq.gz "$i"_L004_R1_001.fastq.gz --validateMappings -o ../../../../salmon/counts/"$i"_out
done;

#multiQC
cd /FDC_rnaseq
multiqc .

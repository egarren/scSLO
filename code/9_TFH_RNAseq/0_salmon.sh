
#salmon index
mkdir salmon
cd salmon
wget http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/salmon_sa_index?tag=default
tar zxvf salmon_sa_index.tgz
mv default mm10_index

#salmon count
cd /data
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
    do salmon quant -i ../../../../salmon/mm10_index -l A -r "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz --validateMappings -o ../../../../salmon/counts/"$i"_out
done;

#multiQC
module load multiqc/1.21
multiqc .

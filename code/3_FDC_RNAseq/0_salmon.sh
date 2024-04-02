
#salmon count
cd /FDC_rnaseq
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
    do salmon quant -i ../../../../salmon/mm10_index -l A -r "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz "$i"_L003_R1_001.fastq.gz "$i"_L004_R1_001.fastq.gz --validateMappings -o ../../../../salmon/counts/"$i"_out
done;

#multiQC
cd /FDC_rnaseq
multiqc .

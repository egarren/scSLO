##### Download data ######

#download SRA data from https://www.ncbi.nlm.nih.gov/Traces/study/
# cd /RNAseq_fastq/SRA2
# aws s3 sync s3://sra2/ .
#split into one folder per sample, make sure all fastqs end in _S1_L001_R1_001 etc

#download EBI data
cd /RNAseq_fastq/
mkdir figshare 
cd figshare
wget -O humBRC_sce.rds https://figshare.com/ndownloader/files/37850658
wget -O folBRCsub_sce.rds https://figshare.com/ndownloader/files/37850634 
wget -O tonFRC_sce.rds https://figshare.com/ndownloader/files/37850823
wget -O allBRCint_sce.rds https://figshare.com/ndownloader/files/37636160
wget -O IMMsub_sce.rds https://figshare.com/ndownloader/files/37850646
wget -O humIMM_sce.rds https://figshare.com/ndownloader/files/37850655
wget -O tonAll_sce.rds https://figshare.com/ndownloader/files/37850820
wget -O tonIMM_sce.rds https://figshare.com/ndownloader/files/37850826


#download EBI data
cd /RNAseq_fastq/EBI
lftp ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/
mirror -P10 E-MTAB-8445 E-MTAB-8445
mirror -P10 E-MTAB-7094 E-MTAB-7094
mirror -P10 E-MTAB-10196 E-MTAB-10196
mirror -P10 E-MTAB-10197 E-MTAB-10197
mirror -P10 E-MTAB-10196 E-MTAB-10196
mirror -P10 E-MTAB-10206 E-MTAB-10206
mirror -P10 E-MTAB-8906 E-MTAB-8906
mirror -P10 E-MTAB-7703 E-MTAB-7703
#split into one folder per sample, make sure all fastqs end in _S1_L001_R1_001 etc

#donwload SCP data
cd /RNAseq_fastq/SCP
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1186&auth_code=F4S46sSy&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1423&auth_code=6RxMdoiu&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2169&auth_code=FOrwNWcF&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1422&auth_code=sJ7E9UCP&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt



###### ALIGNMENTS ######

module load star
cd /RNAseq_fastq/STARsolo
# for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
for d in SRR14355140 SRR14355144 SRR14355148 SRR14355152; do
  cd ./$d
  #concatenate files
  echo "Merging $d R1"
  cat *R1*.fastq.gz > "$d"_R1.fastq.gz
  echo "Merging $d R2"
  cat *R2*.fastq.gz > "$d"_R2.fastq.gz
  echo "Aligning $d"
  STAR --genomeDir /ref_genomes/GRCh38/ref \
  --readFilesIn "$d"_R2.fastq.gz  "$d"_R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBstart 7 \
  --soloCBlen 6  \
  --soloUMIstart 1  \
  --soloUMIlen 6  \
  --soloBarcodeMate 2 \
  --soloCBwhitelist None \
  --clip5pNbases 0 36
  cd ..
done;


#### Cell ranger pipeline (for 10X experiments)
module load cellranger

#process any bam files
cd /RNAseq_fastq/SRA/SRR6252306
cellranger bamtofastq --nthreads=8 ./B5_possorted_genome_bam.bam.1 /RNAseq_fastq/SRA/SRR7227541
cellranger bamtofastq --nthreads=8 ./B6_possorted_genome_bam.bam.1 /RNAseq_fastq/SRA/SRR7227542
cellranger bamtofastq --nthreads=8 ./L1700566_possorted_genome_bam.bam.1 /RNAseq_fastq/SRA/SRR6252306
cellranger bamtofastq --nthreads=8 ./L1700567_possorted_genome_bam.bam.1 /RNAseq_fastq/SRA/SRR6252307

#count matrices
cd /RNAseq_cellranger/FDC

sed -i -e 's/\r$//' sra_mouse.txt
for d in $(cat sra_mouse.txt); do
     cellranger count --id=$d --nosecondary --transcriptome=/ref_genomes/mm39 \
--fastqs=/RNAseq_fastq/SRA/$d
done

sed -i -e 's/\r$//' sra_human.txt
for d in $(cat sra_human.txt); do
     cellranger count --id=$d --nosecondary --transcriptome=/ref_genomes/refdata-gex-GRCh38-2020-A \
--fastqs=/RNAseq_fastq/SRA/$d
done

#if sample is poor quality, must specifcy chemistry
for d in SRR8387860; do
cellranger count --id=$d --nosecondary --transcriptome=/ref_genomes/refdata-gex-GRCh38-2020-A \
--fastqs=/RNAseq_fastq/SRA/$d --chemistry SC3Pv2
done




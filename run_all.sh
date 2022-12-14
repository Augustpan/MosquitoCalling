# TUNE THE FOLLOWING PARAMETERS BEFORE RUNNING
export REF_BT2=./db/BOLD_Culicidae_99_sp_clean
export COI_DB=./db/Cyt_C_Oxase_1_IPR000883.dmnd
export BLASTN_DB=/shilab2/home/public/database/nt/nt
export DIAMOND_EXEC=/usr/local/bin/diamond
export DATA_DIR=/shilab2/home/panyuanfei/mos_rRNA/norRNA_reads
export READ1_SUFFIX=".norRNA.fq.1.gz"
export READ2_SUFFIX=".norRNA.fq.2.gz"
export COI_READS_OUTPUT_DIR="coi_reads"
export COI_CONTIG_OUTPUT_DIR="contigs"
export SNP_OUTPUT_DIR="snp_calling"
export NP_BT2=10
export NP_MEGAHIT=10
export NP_BLASTN=40

# DO NOT CHANGE THE FOLLOWING CODES
mkdir $COI_READS_OUTPUT_DIR
mkdir $COI_CONTIG_OUTPUT_DIR
mkdir $SNP_OUTPUT_DIR

cat manifest.txt \
    | grep -v "^$" \
    | parallel -j4 sh get_COI_reads.sh {} > bowtie2_ouput.log 2>&1

cat manifest.txt \
    | grep -v "^$" \
    | parallel -j4 sh reassemble_COI.sh {} > megahit_ouput.log 2>&1

# skip snp calling
#cat manifest.txt \
#    | grep -v "^$" \
#    | parallel -j4 sh call_snp.sh {} > snpcall_ouput.log 2>&1

#ls $SNP_OUTPUT_DIR/*.vcf \
#    | parallel -j40 python stat_snp.py {} > stat_snp.txt

cat manifest.txt \
    | grep -v "^$" \
    | parallel -j40 seqkit replace -p "^" -r "{}_" "$COI_CONTIG_OUTPUT_DIR/{}.megahit.fa" \
    > all_renamed_contigs.fa

$DIAMOND_EXEC blastx \
    -d $COI_DB \
    -q all_renamed_contigs.fa \
    -o all_renamed_contigs.blastx \
    -e 1e-10 \
    -k 1 \
    -p 40

cut -f 1 all_renamed_contigs.blastx \
    | sort \
    | uniq \
    > coi_contigs.list

seqtk subseq all_renamed_contigs.fa coi_contigs.list > coi_contigs.fa

blastn \
    -query all_renamed_contigs.fa \
    -db $BLASTN_DB \
    -evalue 1e-10 \
    -max_hsps 1 \
    -max_target_seqs 10 \
    -perc_identity 90 \
    -qcov_hsp_perc 30 \
    -num_threads $NP_BLASTN \
    -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle" \
    -out all_contigs_vs_nt.txt

cut -f 13 all_contigs_vs_nt.txt \
    | sort \
    | uniq \
    | taxonkit lineage -c \
    | awk '$2>0' \
    | cut -f 2- \
    | taxonkit reformat -F -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
    | cut -f 1,3- \
    | csvtk join -Ht -L -f "13;1" all_contigs_vs_nt.txt - \
    > all_contigs_vs_nt.parsed.txt

rm all_contigs_vs_nt.txt
mv all_contigs_vs_nt.parsed.txt all_contigs_vs_nt.txt

cut -f 1 -d "_" all_contigs_vs_nt.txt > lib_id
paste lib_id all_contigs_vs_nt.txt \
    | csvtk summary -Ht -g "1,2" -f "3:first,4:first,5:first,12:first,13:first,17:first,18:first,19:first,20:first,21:first,22:first" \
    > all_contigs_tophits.txt
rm lib_id

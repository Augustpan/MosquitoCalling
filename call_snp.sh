LIB_ID=$1

cat "$COI_CONTIG_OUTPUT_DIR/$LIB_ID.megahit.fa" \
    | seqkit sort -l 2> /dev/null \
    | seqkit seq -n \
    | tail -1 \
    > "$LIB_ID.seqname"

seqtk subseq "$COI_CONTIG_OUTPUT_DIR/$LIB_ID.megahit.fa" "$LIB_ID.seqname" \
    > "$LIB_ID.longest.fa"

bowtie2-build "$LIB_ID.longest.fa" "$LIB_ID.longest"

bowtie2 \
    --end-to-end \
    -x "$LIB_ID.longest" \
    -1 "$COI_READS_OUTPUT_DIR/$LIB_ID.1.fq" \
    -2 "$COI_READS_OUTPUT_DIR/$LIB_ID.2.fq" \
    -S "$LIB_ID.sam" \
    --threads $NP_BT2

samtools view -bSF4 $LIB_ID.sam > $LIB_ID.bam
samtools sort $LIB_ID.bam > $SNP_OUTPUT_DIR/$LIB_ID.sorted.bam
samtools index $SNP_OUTPUT_DIR/$LIB_ID.sorted.bam

samtools mpileup -uf "$LIB_ID.longest.fa" $SNP_OUTPUT_DIR/$LIB_ID.sorted.bam 2> /dev/null \
    | bcftools call --ploidy 1 -c \
    > $SNP_OUTPUT_DIR/$LIB_ID.vcf

mv "$LIB_ID.longest.fa" $SNP_OUTPUT_DIR

rm "$LIB_ID.seqname"
rm $LIB_ID.sam $LIB_ID.bam
rm "$LIB_ID".longest.{1,2,3,4,rev.1,rev.2}.bt2
rm "$LIB_ID.longest.fa.fai"

LIB_ID=$1

megahit \
    --num-cpu-threads $NP_MEGAHIT \
    --min-contig-len 300 \
    -1 "$COI_READS_OUTPUT_DIR/$LIB_ID.1.fq" \
    -2 "$COI_READS_OUTPUT_DIR/$LIB_ID.2.fq" \
    -o "$LIB_ID.megahit"

cat "$LIB_ID.megahit/final.contigs.fa" \
    | sed "s/=//g" \
    | sed "s/ /_/g" \
    > "$COI_CONTIG_OUTPUT_DIR/$LIB_ID.megahit.fa"

rm -rf "$LIB_ID.megahit"

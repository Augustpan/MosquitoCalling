LIB_ID=$1

bowtie2 \
    --local \
    -x "$REF_BT2" \
    -1 "$DATA_DIR/${LIB_ID}${READ1_SUFFIX}" \
    -2 "$DATA_DIR/${LIB_ID}${READ2_SUFFIX}" \
    -S "$LIB_ID.sam" \
    --al-conc "$COI_READS_OUTPUT_DIR/$LIB_ID.fq" \
    --threads $NP_BT2

rm "$LIB_ID.sam"

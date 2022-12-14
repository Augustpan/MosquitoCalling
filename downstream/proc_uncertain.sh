# download the following files:
#   all_renamed_contigs.fa
#   all_contigs_tophits.txt
#   coi_contigs.list
#   stat_snp.txt

#Rscript species_calling.R

diamond blastx \
    -d ../db/Cyt_C_Oxase_1_IPR000883.dmnd \
    -q intermediate/all_renamed_contigs.fa \
    -e 1e-10 \
    -k 1 \
    -p 4 \
    -f 6 qseqid qstrand \
    -o intermediate/all_contigs.strand

python convert_strand.py

# extract all uncertain contigs
seqkit grep -r -n \
    -f intermediate/uncertain_libs.txt \
    intermediate/all_coi_contigs_pos.fa \
    > intermediate/uncertain_libs_contigs.fa

# extract mosquito contigs only
seqtk subseq \
    intermediate/uncertain_libs_contigs.fa \
    intermediate/culicidae_contigs.txt \
    > intermediate/coi_mos.fa

# multiple sequences alignment
hmmalign \
    --trim \
    -o intermediate/uncertain_libs_contigs.aln \
    ../db/BOLD_Culicidae_99_sp_clean.hmm \
    intermediate/coi_mos.fa

# calculate pairwise identity and do clustering
python calc_pident_matrix.py \
    intermediate/uncertain_libs_contigs.aln \
    intermediate/uncertain_libs_contigs

cut -f 1 -d "_" \
    intermediate/uncertain_libs_contigs.cluster \
    > intermediate/lib_id

paste \
    intermediate/lib_id \
    intermediate/uncertain_libs_contigs.cluster \
    > intermediate/uncertain_libs_contigs.cluster.edited

cat intermediate/uncertain_libs_contigs.cluster.edited \
    | cut -f 1,3 \
    | sort \
    | uniq \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | sort \
    > intermediate/num_98cluster_per_lib.txt

cat intermediate/num_98cluster_per_lib.txt \
    | sed -r 's/ {2,}//g' \
    | grep -v "^1 " \
    | cut -f2 -d" " \
    > results/contamination_identified_from_uncertain_libs.txt

cat intermediate/num_98cluster_per_lib.txt \
    | sed -r 's/ {2,}//g' \
    | grep "^1 " \
    | cut -f2 -d" " \
    > results/noncontamination_identified_from_uncertain_libs.txt

python update_table.py
python extract_repseq.py

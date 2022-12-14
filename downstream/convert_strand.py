from Bio import SeqIO

with open("intermediate/all_contigs.strand") as f:
    strand_map = [line.split() for line in f.readlines() if line.strip()]

strand_map = dict(strand_map)
ret = []
for rec in SeqIO.parse("intermediate/all_renamed_contigs.fa", "fasta"):
    if rec.id in strand_map:
        if strand_map[rec.id] == "-":
            rec.seq = rec.seq.reverse_complement()
            ret.append(rec)
        else:
            ret.append(rec)

SeqIO.write(ret, "intermediate/all_coi_contigs_pos.fa", "fasta")

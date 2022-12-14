from Bio import SeqIO
import re

with open("summary_table.csv") as f:
    lines = f.readlines()[1:]

single_species_libs = []
mixed_species_libs = []
for line in lines:
    line = line.strip()
    if line:
        sp = line.split(",")
        lib_id = sp[0]
        is_contaminated = sp[3]
        species = sp[2]
        if is_contaminated != "yes":
            single_species_libs.append((lib_id, species))
        else:
            mixed_species_libs.append((lib_id, species))

with open("intermediate/coi_contigs.list") as f:
    coi_contigs = [line.strip() for line in f.readlines() if line.strip()]

seq_dict = {}
for rec in SeqIO.parse("intermediate/all_coi_contigs_pos.fa", "fasta"):
    if rec.id in coi_contigs:
        lib_id = re.findall("(.+?)_k", rec.id).pop()
        if lib_id not in seq_dict:
            seq_dict[lib_id] = []
        seq_dict[lib_id].append(rec)

ret = []
for lib_id, species in single_species_libs:
    longest = sorted(seq_dict[lib_id], key=lambda x: len(x.seq), reverse=True).pop()
    ret.append((lib_id, species, longest))

with open("intermediate/uncertain_libs.txt") as f:
    uncertain_libs = [line.strip() for line in f.readlines() if line.strip()]

with open("intermediate/uncertain_libs_contigs.cluster") as f:
    lines = f.readlines()

cluster_dict = {}
for line in lines:
    line = line.strip()
    if line:
        sp = line.split("\t")
        contig = sp[0]
        cluster = sp[1]
        lib_id = re.findall("(.+?)_k", contig).pop()
        length = int(re.findall("_len(\d+)$", contig).pop())
        if lib_id not in cluster_dict:
            cluster_dict[lib_id] = {}
        if cluster not in cluster_dict[lib_id]:
            cluster_dict[lib_id][cluster] = []
        cluster_dict[lib_id][cluster].append((contig, length))

with open("intermediate/all_contigs_tophits.txt") as f:
    lines = f.readlines()
tophits_dict = {}
for line in lines:
    line = line.strip()
    if line:
        sp = line.split("\t")
        if len(sp) == 13 and sp[10] == "Culicidae":
            lib_id = sp[0]
            contig = sp[1]
            length = int(re.findall("_len(\d+)$", contig).pop())
            species = sp[12]
            if lib_id not in tophits_dict:
                tophits_dict[lib_id] = {}
            if cluster not in tophits_dict[lib_id]:
                tophits_dict[lib_id][species] = []
            tophits_dict[lib_id][species].append((contig, length))

mixed_ret = []
for lib_id, species in mixed_species_libs:
    if lib_id in uncertain_libs:
        for cluster in cluster_dict[lib_id]:
            longest_name, _ = sorted(cluster_dict[lib_id][cluster], key=lambda x: x[1], reverse=True).pop()
            for rec in seq_dict[lib_id]:
                if rec.id == longest_name:
                    mixed_ret.append((lib_id, species, rec))
                    break
    else:
        for spe in tophits_dict[lib_id]:
            longest_name, _ = sorted(tophits_dict[lib_id][spe], key=lambda x: x[1], reverse=True).pop()
            for rec in seq_dict[lib_id]:
                if rec.id == longest_name:
                    mixed_ret.append((lib_id, spe, rec))
                    break

ret_seqs = []
for lib_id, species, rec in ret:
    rec.description = f"[lib_id: {lib_id}] [species: {species}]"
    ret_seqs.append(rec)

ret_seqs_mixed = []
for lib_id, species, rec in mixed_ret:
    rec.description = f"[lib_id: {lib_id}] [species: {species}]"
    ret_seqs_mixed.append(rec)

SeqIO.write(ret_seqs, "results/repseqs_single.fa", "fasta")
SeqIO.write(ret_seqs_mixed, "results/repseqs_mixed.fa", "fasta")

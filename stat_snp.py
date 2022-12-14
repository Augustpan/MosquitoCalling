import pandas as pd
import re
import sys
import pathlib

df = pd.read_csv(sys.argv[1], sep="\t", comment="#", header=None)

ret = map(lambda x: re.findall("DP4=(\d+),(\d+),(\d+),(\d+);", x)[0], df[7])
ref_v_var = list(map(lambda x: (int(x[0])+int(x[1]), int(x[2])+int(x[3])), ret))

seq_len = int(re.findall("_len(\d+)$", df[0][0]).pop())
snp_site_count = 0
mean_freq = 0
mean_depth = 0
for ref, var in ref_v_var:
    if var == 0 or ref == 0:
        continue
    tot = sum([ref+var])
    var_freq = min([ref, var]) / tot
    var_depth = min([ref, var])
    if var_depth > 5 and var_freq > 0.05:
        snp_site_count += 1
        mean_freq += var_freq
        mean_depth += var_depth
if snp_site_count:
    mean_freq /= snp_site_count
    mean_depth /= snp_site_count
    snp_site_count /= seq_len

print(f"{pathlib.Path(sys.argv[1]).stem}\t{snp_site_count* 1000}\t{mean_freq}\t{mean_depth}")

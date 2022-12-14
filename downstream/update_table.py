with open("intermediate/all.csv") as f:
    lines = f.readlines()

with open("results/contamination_identified_from_uncertain_libs.txt") as f:
    cont = [line.strip() for line in f.readlines() if line.strip()]

new_lines = []


for line in lines:
    line = line.strip()
    if line:
        sp = line.split(",")
        if sp[3] == "yes":
            if sp[2] == "Culex pipiens/Culex quinquefasciatus":
                sp[3] = "indistinguishable"
        elif sp[3] == "unsure":
            for cc in cont:
                if sp[0] == cc:
                    sp[3] = "yes"
                    break
            if sp[3] == "unsure":
                sp[3] = "no_but_species_uncertain"
        new_lines.append(",".join(sp))
with open("results/summary_table.csv", "w") as f:
    f.write("\n".join(new_lines) + "\n")

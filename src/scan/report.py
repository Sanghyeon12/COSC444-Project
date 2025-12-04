# report.py – Sadman (Person 3)

def write_matches(out_path, results):
    f = open(out_path, "w")
    for seq_id, hits in results.items():
        for h in hits:
            line = seq_id + "\t" + h[0] + "\t" + str(h[1]) + "\t" + str(h[2]) + "\t" + h[3] + "\n"
            f.write(line)
    f.close()

def write_summary(out_path, results):
    f = open(out_path, "w")
    for seq_id, hits in results.items():
        f.write(seq_id + "\t" + str(len(hits)) + "\n")
    f.close()

def write_baseline(out_path, results):
    f = open(out_path, "w")
    for seq_id, hits in results.items():
        f.write(seq_id + "\t" + str(len(hits)) + "\n")
    f.close()

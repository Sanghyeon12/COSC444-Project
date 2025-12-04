import re

def reverse_comp(seq):
    comp = {"A":"T","T":"A","C":"G","G":"C",
            "a":"t","t":"a","c":"g","g":"c"}
    out = []
    for s in seq:
        if s in comp:
            out.append(comp[s])
        else:
            out.append(s)
    out.reverse()
    return "".join(out)

def find_matches(seq, patterns, allow_overlap):
    hits = []
    for name, pat in patterns.items():
        if allow_overlap:
            for m in re.finditer(pat, seq):
                hits.append((name, m.start(), m.end(), m.group(), "+"))
        else:
            pos = 0
            while True:
                m = pat.search(seq, pos)
                if not m:
                    break
                hits.append((name, m.start(), m.end(), m.group(), "+"))
                pos = m.end()
    return hits

def scan_sequence(seq, patterns, allow_overlap, revcomp):
    all_hits = []
    forward = find_matches(seq, patterns, allow_overlap)
    all_hits.extend(forward)

    if revcomp:
        rc = reverse_comp(seq)
        rc_hits = find_matches(rc, patterns, allow_overlap)
        for h in rc_hits:
            all_hits.append((h[0], h[1], h[2], h[3], "-"))

    return all_hits

def scan_fasta(records, patterns, allow_overlap, revcomp):
    result = {}
    for name, seq in records.items():
        result[name] = scan_sequence(seq, patterns, allow_overlap, revcomp)
    return result

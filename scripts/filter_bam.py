import subprocess, sys, re
src, dst = sys.argv[1], sys.argv[2]
chroms = {l.split("\t")[0] for l in open("/home/dbs/dna_data/raw/GRCh38_no_alt.fa.fai")}
p = subprocess.Popen(["samtools", "view", "-h", src], stdout=subprocess.PIPE, text=True)
q = subprocess.Popen(["samtools", "view", "-b", "-o", dst, "-"], stdin=subprocess.PIPE, text=True)
n = d = 0
for l in p.stdout:
    if l[0] == "@":
        q.stdin.write(l); continue
    f = l.split("\t")
    n += 1
    if f[6] not in ("=", "*"):
        d += 1; continue
    tags = []
    for t in f[11:]:
        if t.startswith("SA:Z:"):
            segs = [s for s in t[5:].rstrip("\n").split(";") if s and s.split(",")[0] in chroms]
            if not segs: continue
            t = "SA:Z:" + ";".join(segs) + ";"
        tags.append(t.rstrip("\n"))
    q.stdin.write("\t".join(f[:11] + tags) + "\n")
q.stdin.close(); q.wait()
print(src, n, "reads;", d, "dropped (mate off-reference)")

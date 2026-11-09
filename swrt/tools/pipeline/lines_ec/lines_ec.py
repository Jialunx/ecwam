"""Searches along the line of 4 times the cost of DIA on ecWAM's basis (../cost_ec/time_sets_ec.json: SNONLIN_SWRT and SNONLIN of
bundle/build_bit on the 15,000 synthetic training spectra; t = c0 + a K + b L), exactly as keydoc/swrt/candidates/ww3cost/scripts/
lines_ww3.py did on WAVEWATCH III's basis: SWRT's search at full effort (candidates/K36/scripts/search_k36.one_L of this grid: greedy
under the offset limit, swaps, 6 random starts, 30 kicks; seeds 11-13) for every K on the line with 6 <= L <= 2K,
L = floor((c DIA - c0 - a K) / b), on the summed synthetic Gram of this grid (results/gram_synthetic_s*.npz, 10 draws). A search depends
only on (K, L, seed), so full-effort results already on this grid (candidates/K36/results/search_k36.json, K = 36) are reused when
their (K, L, seed) is on the line. Appends line_runs_ec.jsonl (one record per line, K and seed; resumable: records already in the
file are skipped). Tools copy: given a list of K, only those K of the line are searched (the tool needs K = 31 only).
usage: OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 lines_ec.py 4 [workers] [K,K,...]"""
import os, sys, json, time, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "candidates", "K36", "scripts"))
T = json.load(open(os.path.join(ROOT, "cost_ec", "time_sets_ec.json"))); FT, DIA = T["fit"], T["dia_us"]
cost_of = lambda K, L: (FT["c0"] + FT["per_K"] * K + FT["per_L"] * L) / DIA
L_of = lambda c, K: int(np.floor((c * DIA - FT["c0"] - FT["per_K"] * K) / FT["per_L"] + 1e-9))
OUT = os.path.join(HERE, "line_runs_ec.jsonl")


def known():
    H = {}
    K36 = os.path.join(ROOT, "candidates", "K36", "results", "search_k36.json")       # tools copy: absent here
    for k, r in (json.load(open(K36)) if os.path.exists(K36) else {}).items():
        if k != "note": H[(36, r["L"], r["seed"])] = dict(units=r["units"], e2=r["e2"], loci=r["loci"], src="search_k36.json")
    if os.path.exists(OUT):
        for r in map(json.loads, open(OUT)):
            if r.get("e2") is not None: H.setdefault((r["K"], r["L"], r["seed"]), dict(units=r["units"], e2=r["e2"], loci=r["loci"], src="line_runs_ec.jsonl"))
    return H


def job(a):
    c, K, L, seed = a; t0 = time.time(); import search_k36 as sk
    try: u, e = sk.one_L(L, np.random.default_rng(seed), K=K)
    except ValueError as x:                                                     # as lines_ww3.py: a set too small for the kicks, or no feasible set
        return dict(line=c, K=K, L=L, seed=seed, starts=6, kicks=30, units=[], e2=None, loci=None, cost=None, error=repr(x)[:200], seconds=time.time() - t0, src="run")
    lo = int(sk.lb.loci(u))
    return dict(line=c, K=K, L=L, seed=seed, starts=6, kicks=30, units=[int(x) for x in u], e2=float(e), loci=lo, cost=cost_of(K, lo), seconds=time.time() - t0, src="run")


if __name__ == "__main__":
    lines = [float(x) for x in sys.argv[1].split(",")]; nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8; H = known()
    KS = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else range(6, 120)
    have = {(r["line"], r["K"], r["seed"]) for r in map(json.loads, open(OUT))} if os.path.exists(OUT) else set()
    jobs, reuse = [], []
    for c in lines:
        for K in KS:
            L = L_of(c, K)
            if not (6 <= L <= 2 * K): continue
            for s in (11, 12, 13):
                if (c, K, s) in have: continue
                if (K, L, s) in H: h = H[(K, L, s)]; reuse.append(dict(line=c, K=K, L=L, seed=s, starts=6, kicks=30, units=h["units"], e2=h["e2"], loci=h["loci"], cost=cost_of(K, h["loci"]), seconds=0.0, src=h["src"]))
                else: jobs.append((c, K, L, s))
    print("lines", lines, ":", {c: [(K, L_of(c, K)) for K in KS if 6 <= L_of(c, K) <= 2 * K] for c in lines}, flush=True)
    print(len(reuse), "reused,", len(jobs), "to run", flush=True)
    with open(OUT, "a") as fh:
        for r in reuse: fh.write(json.dumps(r) + "\n")
    jobs.sort(key=lambda j: -j[1] * j[2]); t0 = time.time()
    with Pool(nw) as p, open(OUT, "a") as fh:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            fh.write(json.dumps(r) + "\n"); fh.flush()
            print(f"{i + 1}/{len(jobs)} K{r['K']} L{r['L']} seed {r['seed']}: " + (f"e2 {r['e2']:.5f} uses {r['loci']}" if r["e2"] is not None else r["error"])
                  + f" ({r['seconds']:.0f} s; {time.time() - t0:.0f} s)", flush=True)
    R = [r for r in map(json.loads, open(OUT)) if r.get("e2") is not None]
    for c in lines:
        Rc = [r for r in R if r["line"] == c]; byK = {}
        for r in Rc:
            if r["K"] not in byK or r["e2"] < byK[r["K"]]["e2"]: byK[r["K"]] = r
        for K in sorted(byK): print(f"  {c:g}x K {K:2d} L {byK[K]['L']:2d} (uses {byK[K]['loci']:2d}) e2 {byK[K]['e2']:.5f} cost {byK[K]['cost']:.2f}x seed {byK[K]['seed']}")
        b = min(Rc, key=lambda r: r["e2"]); print(f"best on {c:g}x: K{b['K']}/{b['loci']} (limit {b['L']}, seed {b['seed']}) e2 {b['e2']:.5f} cost {b['cost']:.2f}x", flush=True)

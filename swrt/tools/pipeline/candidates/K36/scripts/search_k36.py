"""The search that keeps SWRT's configurations since 2026-10-20 (late): K = 36 configurations whose rows use at most L = 18 offsets
(di, dj), the set with the smallest ebar_2 on the summed synthetic Gram at a kernel cost of 2.5 times DIA (user: "use 36/18 for my new
best swrt"). one_L is the search of candidates/N30/scripts/locibudget.py (greedy under the offset limit, swaps, random feasible starts,
kicks), run with numpy's default_rng(seed) for one limit at a time: SWRT is one_L(18, default_rng(12), K = 36). Also run here, for
section 3b: the limits 14 and 16 (seed 12) and no limit (L = 999, seeds 11, 12, 13).
usage: OMP_NUM_THREADS=1 python3 search_k36.py  ->  ../results/search_k36.json"""
import os, sys, json, time, numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "candidates", "N30", "scripts")); import locibudget as lb


def one_L(L, rng, K=36, starts=6, kicks=30):
    best = lb.swap(lb.greedy(K, L), L); eb = lb.e2(best)
    for _ in range(starts):
        Sx = []
        for q in rng.permutation(lb.N):
            if len(Sx) == K: break
            if lb.loci(Sx + [int(q)]) <= L: Sx.append(int(q))
        if len(Sx) < K: continue
        Sx = lb.swap(Sx, L)
        if lb.e2(Sx) < eb - 1e-12: best, eb = Sx, lb.e2(Sx)
    for _ in range(kicks):
        Sx = list(best)
        for pos in rng.choice(K, int(rng.integers(3, 7)), replace=False):
            c = [q for q in range(lb.N) if q not in Sx and lb.loci(Sx[:pos] + Sx[pos + 1:] + [q]) <= L]
            if c: Sx[pos] = int(rng.choice(c))
        Sx = lb.swap(Sx, L)
        if lb.e2(Sx) < eb - 1e-12: best, eb = Sx, lb.e2(Sx)
    return sorted(int(u) for u in best), float(eb)


def job(a):
    L, seed = a; t0 = time.time(); u, e = one_L(L, np.random.default_rng(seed))
    return dict(L=L, seed=seed, units=u, e2=e, train=100 * np.sqrt(e), loci=lb.loci(u), seconds=time.time() - t0)


if __name__ == "__main__":
    jobs = [(18, 12), (14, 12), (16, 12), (999, 11), (999, 12), (999, 13)]
    with Pool(6) as p: res = p.map(job, jobs)
    out = {f"{r['L']}_{r['seed']}": r for r in res}
    opf = os.path.join(ROOT, "operators", "swrt.npz")              # ecWAM copy: written from 18_12 by make_swrt_op_K36L18.py afterwards
    sw = sorted(int(u) for u in np.load(opf)["units"]) if os.path.exists(opf) else None
    out["note"] = "K = 36; one_L(L, default_rng(seed)); SWRT (operators/swrt.npz) = 18_12: " + str(out["18_12"]["units"] == sw if sw else "operator not yet written")
    json.dump(out, open(os.path.join(HERE, "..", "results", "search_k36.json"), "w"), indent=1)
    for k, r in out.items():
        if k != "note": print(k, f"e2 {r['e2']:.5f} offsets {r['loci']} {r['seconds']:.0f} s")
    print(out["note"])

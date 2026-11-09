"""SWRT sets for a speed budget (the user, 2026-11-08: "if I want 1.5 x the speed, give me the most accurate and most stable one under
that limit"; K and the set chosen from training data only). SWRT's objective and setting: e2 on the summed synthetic Gram of the
practice pool (seeds 11, 1-9), least-squares gains (positive), c = 1; the search of search_syn.py with e2 only (greedy under the offset
limit, swaps, 6 random starts, 30 kicks), three search seeds (11-13), for every K = 2-12 and every offset limit L with the modelled
S_nl cost (cost_ec/time_sets_ec.json: t = c0 + a K + b L, over DIA's time) at most the budget. The best set per (K, L limit) is kept.
usage: python3 search_budget.py <grid root> <budget, x DIA> [workers]  ->  <grid root>/synfast/budget_<budget>.json (never overwritten)"""
import os, sys, json, time, numpy as np
from multiprocessing import Pool
ROOT, BUD = os.path.abspath(sys.argv[1]), float(sys.argv[2]); NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
OUT = os.path.join(ROOT, "synfast", f"budget_{BUD:g}.json")
FT = json.load(open(os.path.join(ROOT, "cost_ec", "time_sets_ec.json"))); C0, AK, BL, DIA = FT["fit"]["c0"], FT["fit"]["per_K"], FT["fit"]["per_L"], FT["dia_us"]
cost = lambda K, L: (C0 + AK * K + BL * L) / DIA
T = dict(np.load(f"{ROOT}/data/table.npz")); NU = int(T["unit"].max()) + 1
pair = [np.flatnonzero(T["unit"] == u) for u in range(NU)]
LOC = [frozenset((int(T["di"][r]), int(T["dj"][r])) for r in pair[u]) for u in range(NU)]
loci = lambda S: len(frozenset().union(*[LOC[u] for u in S])) if S else 0
G = 0; b = 0; yy = 0.0
for s in (11, 1, 2, 3, 4, 5, 6, 7, 8, 9):
    z = np.load(f"{ROOT}/results/gram_synthetic_s{s}.npz"); G = G + z["G"]; b = b + z["b"]; yy += float(z["yy"])
Gs, bs_ = G / yy, b / yy


def fit(S): g = np.linalg.solve(Gs[np.ix_(S, S)], bs_[S]); return g, float(1 - bs_[S] @ g)


def J(S):
    if not S: return 1.0
    try: g, j = fit(S)
    except np.linalg.LinAlgError: return np.inf
    return j if (g > 0).all() else np.inf


def greedy(K, L):
    S = []
    while len(S) < K:
        c = [q for q in range(NU) if q not in S and loci(S + [q]) <= L]
        if not c: return None
        v = [J(S + [q]) for q in c]
        if not np.isfinite(min(v)): return None
        S.append(c[int(np.argmin(v))])
    return S


def swap(S, L):
    S = list(S)
    while True:
        best = J(S)
        for pos in range(len(S)):
            rest = S[:pos] + S[pos + 1:]; c = [q for q in range(NU) if q not in rest and loci(rest + [q]) <= L]
            v = [J(rest + [q]) for q in c]; q = c[int(np.argmin(v))]
            if min(v) < J(S) - 1e-12: S = rest + [q]
        if J(S) >= best - 1e-12: return S


def search(K, L, rng, starts=6, kicks=30):
    best = greedy(K, L); best = swap(best, L) if best else None; eb = J(best) if best else np.inf
    for _ in range(starts):
        Sx = []
        for q in rng.permutation(NU):
            if len(Sx) == K: break
            if loci(Sx + [int(q)]) <= L and np.isfinite(J(Sx + [int(q)])): Sx.append(int(q))
        if len(Sx) < K: continue
        Sx = swap(Sx, L)
        if J(Sx) < eb - 1e-12: best, eb = Sx, J(Sx)
    if best is None: return None, np.inf
    for _ in range(kicks):
        Sx = list(best)
        for pos in rng.choice(K, int(rng.integers(1, max(2, K // 2) + 1)), replace=False):
            c = [q for q in range(NU) if q not in Sx and loci(Sx[:pos] + Sx[pos + 1:] + [q]) <= L]
            if c: Sx[pos] = int(rng.choice(c))
        if not np.isfinite(J(Sx)): continue
        Sx = swap(Sx, L)
        if J(Sx) < eb - 1e-12: best, eb = Sx, J(Sx)
    return sorted(int(u) for u in best), float(eb)


def job(a):
    K, L, seed = a; t0 = time.time(); S, e = search(K, L, np.random.default_rng(seed))
    return dict(K=K, L_limit=L, seed=seed, units=S, e2=e if S else None, L=loci(S) if S else None, seconds=time.time() - t0)


if __name__ == "__main__":
    if os.path.exists(OUT): sys.exit(f"{OUT} exists")
    jobs = [(K, L, s) for K in range(2, 13) for L in range(1, 2 * K + 1) if cost(K, L) <= BUD + 1e-9 for s in (11, 12, 13)]
    print(f"budget {BUD:g} x DIA: {len(jobs) // 3} (K, L) pairs, {len(jobs)} searches", flush=True); t0 = time.time(); R = []
    with Pool(NW) as p:
        for r in p.imap_unordered(job, sorted(jobs, key=lambda j: -j[0] * j[1])):
            R.append(r)
    best = {}
    for r in R:
        if r["units"] is None: continue
        k = (r["K"], r["L_limit"])
        if k not in best or (r["e2"], r["seed"]) < (best[k]["e2"], best[k]["seed"]): best[k] = r
    sets = {}
    for (K, Ll), r in sorted(best.items()):
        g, _ = fit(r["units"]); sets[f"K{K}Llim{Ll}"] = dict(r, g=[float(x) for x in g], cost=cost(K, r["L"]))
    json.dump(dict(budget=BUD, cost_fit=FT["fit"], dia_us=DIA, sets=sets, runs=R, note=__doc__.split("usage")[0]), open(OUT, "w"), indent=1)
    for k, r in sorted(sets.items(), key=lambda kv: kv[1]["e2"]):
        print(f"{k:12s} K{r['K']}/L{r['L']} cost {r['cost']:.2f}x e2 {r['e2']:.4f} seed {r['seed']}", flush=True)
    print(f"{time.time() - t0:.0f} s", flush=True)

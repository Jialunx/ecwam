"""Cost-aware selection for the locus-grouped kernel: K units whose rows use at most L distinct (di, dj) loci.
Cost model of swrt_group_d (-O3, bench_group.py): about 0.52 us per row and 1.07 us per locus per spectrum, so with K fixed
fewer loci are faster. Search: greedy under the budget, then swap passes and kicks that keep the budget (as multistart.py).
Objective and data as build_swrt.py (summed synthetic Gram, 15,000 spectra).
usage: venv/bin/python locibudget.py 30 16,18,20,22 [starts] [kicks]"""
import os, sys, json, time, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))   # tools copy: the work folder
G = 0; b = 0; yy = 0.0
for s in (11, 1, 2, 3, 4, 5, 6, 7, 8, 9):
    z = np.load(f"{ROOT}/results/gram_synthetic_s{s}.npz"); G = G + z["G"]; b = b + z["b"]; yy += float(z["yy"])
G, b = G / yy, b / yy; N = len(b)
T = dict(np.load(f"{ROOT}/data/table.npz")); pair = [np.flatnonzero(T["unit"] == u) for u in range(N)]
LOC = [frozenset((int(T["di"][r]), int(T["dj"][r])) for r in pair[u]) for u in range(N)]
loci = lambda S: len(frozenset().union(*[LOC[u] for u in S])) if S else 0
e2 = lambda S: float(1 - b[S] @ np.linalg.solve(G[np.ix_(S, S)], b[S]))
err = lambda v: 100 * np.sqrt(max(v, 0.0))


def greedy(K, L):
    S = []
    while len(S) < K:
        c = [q for q in range(N) if q not in S and loci(S + [q]) <= L]; S.append(min(c, key=lambda q: e2(S + [q])))
    return S


def swap(S, L):
    S = list(S)
    while True:
        best = e2(S)
        for pos in range(len(S)):
            rest = S[:pos] + S[pos + 1:]; c = [q for q in range(N) if q not in rest and loci(rest + [q]) <= L]
            q = min(c, key=lambda q: e2(rest + [q]))
            if e2(rest + [q]) < e2(S) - 1e-12: S = rest + [q]
        if e2(S) >= best - 1e-12: return S


if __name__ == "__main__":
    K = int(sys.argv[1]); Ls = [int(x) for x in sys.argv[2].split(",")]; starts = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    kicks = int(sys.argv[4]) if len(sys.argv) > 4 else 30; rng = np.random.default_rng(11); out = {}
    for L in Ls:
        t0 = time.time(); best = swap(greedy(K, L), L); eb = e2(best)
        for _ in range(starts):                                           # random feasible starts
            S = []
            for q in rng.permutation(N):
                if len(S) == K: break
                if loci(S + [int(q)]) <= L: S.append(int(q))
            if len(S) < K: continue
            S = swap(S, L)
            if e2(S) < eb - 1e-12: best, eb = S, e2(S)
        for _ in range(kicks):
            S = list(best)
            for pos in rng.choice(K, int(rng.integers(3, 7)), replace=False):
                c = [q for q in range(N) if q not in S and loci(S[:pos] + S[pos + 1:] + [q]) <= L]
                if c: S[pos] = int(rng.choice(c))
            S = swap(S, L)
            if e2(S) < eb - 1e-12: best, eb = S, e2(S)
        out[L] = dict(units=sorted(int(u) for u in best), loci=loci(best), train=err(eb), seconds=time.time() - t0)
        print(f"K={K} L<={L}: loci {loci(best)}, train error {err(eb):.3f}%  ({time.time() - t0:.0f} s)", flush=True)
    json.dump(out, open(f"locibudget_K{K}.json", "w"), indent=1)

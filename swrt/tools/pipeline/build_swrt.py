"""Build the sWRT operator: keep K mirror-pair units of the table, one gain each. No WRT call is needed.

Rows and gains: on a set of practice (probe) spectra, find K units and gains g so that sum_u g_u Psi_u(E) reproduces the
full-table sum sum_all Psi(E). The full-table sum is WRT (to single-precision round-off, check_table.py), so the gains carry
WRT's absolute units and the operator is used with c = 1. Each probe spectrum is scaled by 1 / |full sum| (Snl ~ E^3,
otherwise a few large seas dominate). Search: greedy from empty (each step tries every unit not yet kept, adds the one that
lowers the residual sum of squares most, refits all gains by ordinary least squares), then swap passes (replace one kept
unit at a time by the best other unit) until no swap lowers the residual. The residual is closed form from the Gram matrix.
Check only (not applied): c = <P, S_WRT> / <P, P> from one WRT call on one JONSWAP spectrum.

probe sets: synthetic (default; make_synthetic.py, formulas + seeds, no model data; --seeds pools several draws, 1500 each),
ww3 (data/probe_1200.npz), both (ww3 + synthetic seed --seed).
The operator used in the write-up is K = 30 from the 10-seed pool (15000 spectra): robust to the draw (robustness.py).
usage: python3 build_swrt.py                          ->  operators/swrt_K{10,20,30,50}.npz, results/build.json  (pool, 10 seeds)
       python3 build_swrt.py --probes ww3 --tag ww3 --ks 30   ->  operators/sensitivity/swrt_ww3_K30.npz, results/build_ww3.json
       python3 build_swrt.py --seeds 3 --tag syn3 --ks 20,30
       python3 build_swrt.py --ks 1,5 --json build_smallK.json   ->  operators/swrt_K{1,5}.npz, results/build_smallK.json (build.json kept)"""
import os, json, time, argparse, numpy as np
import swrt
from wrt_ref import wrt, FREQ, DIRS

HERE = os.path.dirname(os.path.abspath(__file__)); KMAX, BATCH = 50, 16
T = swrt.load_table(); R = swrt.Rows(T); NU = int(T["unit"].max()) + 1
pair = np.array([np.flatnonzero(T["unit"] == u) for u in range(NU)])            # (NU, 2) row ids of each unit


def unit_responses(E):
    r = R.responses(E); r = r.reshape(r.shape[0], r.shape[1], -1)
    return r[:, pair[:, 0]] + r[:, pair[:, 1]]                                   # (B, NU, 960)


def jonswap(fp=0.14, hs=2.0, gamma=3.3, s2=20.0, th0=187.5):
    """JONSWAP in f, sigma_a/b = 0.07/0.09, times D(theta) ~ cos^s2((theta - th0)/2), scaled to Hs"""
    f = FREQ; df = np.gradient(f); sg = np.where(f <= fp, 0.07, 0.09)
    S = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(-((f - fp) ** 2) / (2 * (sg * fp) ** 2))
    D = np.abs(np.cos(np.deg2rad(DIRS - th0) / 2)) ** s2; D /= D.sum() * 2 * np.pi / len(DIRS)
    E = S[:, None] * D[None]; return E * (hs ** 2 / 16) / (E.sum(1) * 2 * np.pi / len(DIRS) * df).sum()


def probe_set(kind, seed=11):
    syn = lambda: np.load(os.path.join(HERE, "data", f"synthetic_s{seed}.npz"))["E"].astype(float)
    ww3 = lambda: np.load(os.path.join(HERE, "data", "probe_1200.npz"))["E"].astype(float)
    return {"synthetic": syn, "ww3": ww3, "both": lambda: np.concatenate([ww3(), syn()])}[kind]()


def gram(E, cache):
    """normalised Gram matrix of the unit responses and their projections on the full-table sum (cached)"""
    if os.path.exists(cache):
        z = np.load(cache); return z["G"], z["b"], float(z["yy"]), float(z["seconds"])
    G = np.zeros((NU, NU)); b = np.zeros(NU); yy = 0.0; t0 = time.time()
    for s in range(0, len(E), BATCH):
        U = unit_responses(E[s:s + BATCH]); full = U.sum(1); sc = 1 / np.sqrt((full ** 2).sum(1))
        A = (U * sc[:, None, None]).transpose(1, 0, 2).reshape(NU, -1); y = (full * sc[:, None]).ravel()
        G += A @ A.T; b += A @ y; yy += float(y @ y)
    sec = time.time() - t0; np.savez(cache, G=G, b=b, yy=yy, seconds=sec); print(f"Gram over {len(E)} probe spectra, {NU} units: {sec:.0f} s")
    return G, b, yy, sec


def select(G, b, yy, Ks):
    fit = lambda Ss: np.linalg.solve(G[np.ix_(Ss, Ss)], b[Ss]); sse = lambda Ss: yy - b[Ss] @ fit(Ss)   # OLS residual
    S, path = [], []
    while len(S) < max(Ks):                                         # greedy over all units
        S.append(min((q for q in range(NU) if q not in S), key=lambda q: sse(S + [q]))); path.append(1 - sse(S) / yy)

    def swap(S):
        S = list(S)
        while True:
            best = sse(S)
            for pos in range(len(S)):
                rest = S[:pos] + S[pos + 1:]; q = min((q for q in range(NU) if q not in rest), key=lambda q: sse(rest + [q]))
                if sse(rest + [q]) < sse(S) - 1e-12 * yy: S = rest + [q]
            if sse(S) >= best - 1e-12 * yy: return S
    out = {}
    for K in Ks:
        keep = swap(S[:K]); out[K] = (np.array(keep), fit(keep), 1 - sse(keep) / yy, path[K - 1])
    return out, path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--probes", default="synthetic", choices=("synthetic", "ww3", "both"))
    ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9"); ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--tag", default=""); ap.add_argument("--ks", default="10,20,30,50"); ap.add_argument("--json", default="")
    a = ap.parse_args(); Ks = tuple(int(k) for k in a.ks.split(",")); seeds = [int(x) for x in a.seeds.split(",")]
    if a.probes == "synthetic":                                       # pooled draws: sum of the per-draw Grams
        G, b, yy, gram_s, n = 0, 0, 0.0, 0.0, 0
        for sd in seeds:
            Gi, bi, yi, si = gram(probe_set("synthetic", sd), os.path.join(HERE, "results", f"gram_synthetic_s{sd}.npz"))
            G, b, yy, gram_s, n = G + Gi, b + bi, yy + yi, gram_s + si, n + 1500
        name = "synthetic_s" + "-".join(map(str, seeds))
    else:
        name = a.probes + (f"_s{a.seed}" if a.probes == "both" else ""); E = probe_set(a.probes, a.seed); n = len(E)
        G, b, yy, gram_s = gram(E, os.path.join(HERE, "results", f"gram_{name}.npz"))
    t0 = time.time(); sel, path = select(G, b, yy, Ks); search_s = time.time() - t0
    Ej = jonswap(); Yj = wrt(Ej).ravel(); Uj = unit_responses(Ej)[0]            # the check (one WRT call), not applied
    odir = os.path.join(HERE, "operators") if not a.tag else os.path.join(HERE, "operators", "sensitivity")
    os.makedirs(odir, exist_ok=True)
    out = {"probes": a.probes, "seeds": seeds if a.probes == "synthetic" else None, "n_probe": n, "probe_fit_r2_greedy_path": path,
           "units_total": NU, "rows_total": len(T["di"]), "gram_seconds": gram_s, "search_seconds": search_s,
           "c_check_jonswap": dict(fp=0.14, hs=2.0, gamma=3.3, s2=20)}
    for K, (keep, g, r2, r2g) in sel.items():
        Pj = g @ Uj[keep]; c_chk = float((Pj * Yj).sum() / (Pj * Pj).sum()); rows = pair[keep].ravel()
        fn = f"swrt_K{K}.npz" if not a.tag else f"swrt_{a.tag}_K{K}.npz"
        np.savez(os.path.join(odir, fn), K=K, units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=c_chk,
                 probes=a.probes, seeds=np.array(seeds if a.probes == "synthetic" else [a.seed]), probe_fit_r2=r2, **{k: T[k][rows] for k in swrt.GEOM + ("amem", "node")}, **{k: T[k] for k in swrt.CONST},
                 note="sWRT: rows of WRT's quadrature (2 per unit: row + mirror), gains fitted to the full-table sum on probe spectra, c = 1")
        out[f"K{K}"] = dict(file=os.path.relpath(os.path.join(odir, fn), HERE), units=keep.tolist(), gains=g.tolist(), c_check=c_chk, probe_fit_r2=r2,
                            greedy_only_r2=r2g, neg_gains=int((g < 0).sum()), configs=[(int(T["di"][pair[u, 0]]), int(abs(T["dj"][pair[u, 0]]))) for u in keep])
        print(f"{name} K={K:2d}: probe-fit R2 {r2:.4f} (greedy alone {r2g:.4f}), gains {g.min():.3g}..{g.max():.3g} (negative {(g < 0).sum()}), c check {c_chk:.4f}")
    json.dump(out, open(os.path.join(HERE, "results", a.json or f"build{'_' + a.tag if a.tag else ''}.json"), "w"), indent=1)
    print(f"search {search_s:.1f} s")

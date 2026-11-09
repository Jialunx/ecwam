"""SWRT of this grid by the paper's rule on ecWAM's cost basis (user 2026-11-03): the set with the smallest training ebar_2 among the
searches on the 4x DIA line of ../cost_ec/time_sets_ec.json (line_runs_ec.jsonl, line 4, every K and seeds 11-13), with least-squares
gains on the summed synthetic Gram (10 draws, 15,000 spectra) and c = 1, written as make_swrt_op_K36L18.py writes operators/swrt.npz
(c_check from one WRT call on one JONSWAP spectrum, a check only). Outputs (new files; an existing file is never overwritten):
  swrt_ec4x_K{K}L{L}.npz                               Python operator (swrt.SWRT), L = offsets used
  ../ww3/nl6/swrt_ec4x_K{K}L{L}.ww3op                  operator file for SNONLIN_SWRT (ww3/nl6/export_operator.py format)
  swrt_ec4x.json                                       the choice, the per-(K, seed) results of the line, the offsets and gains
With a K given (the user, 2026-11-08: the same K as the 36 x 36 file, K = 31), the set of that K on the line instead, and
swrt_ec4x_K{K}.json in place of swrt_ec4x.json. Tools copy: with a rank given as well, the set of that rank among the distinct sets of
that K in the same order (rank 0 = the set above), for the integration filter of make_swrt_files.sh when a set blows up.
usage: python3 make_swrt_ec4x.py [K [rank]]"""
import os, sys, json, tempfile, shutil, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [ROOT, os.path.join(ROOT, "ww3", "nl6")]
import swrt, export_operator as eo, wrt_ref, grid
from wrt_ref import wrt
from build_swrt import jonswap, unit_responses, pair, T
XNL = tempfile.mkdtemp(prefix="xnl_", dir=HERE); os.symlink(os.path.join(ROOT, "quad99999.bqf"), os.path.join(XNL, "quad99999.bqf"))
wrt_ref.HERE = XNL                          # Xnl runs (and writes its xnl4v5.* logs) in a scratch folder, not in swrt_op36 itself

R = [r for r in map(json.loads, open(os.path.join(HERE, "line_runs_ec.jsonl"))) if r["line"] == 4.0]
KF = int(sys.argv[1]) if len(sys.argv) > 1 else None; RANK = int(sys.argv[2]) if len(sys.argv) > 2 else 0
ok = [r for r in R if r["e2"] is not None and KF in (None, r["K"])]; seen, dist = set(), []
for r in sorted(ok, key=lambda r: (r["e2"], r["K"], r["seed"])):              # distinct sets in e2 order (seeds may find the same set)
    if tuple(sorted(r["units"])) not in seen: seen.add(tuple(sorted(r["units"]))); dist.append(r)
if RANK >= len(dist): shutil.rmtree(XNL); sys.exit(f"no set of rank {RANK}: {len(dist)} distinct sets of K = {KF}")
best = dist[RANK]
seeds = [11, 1, 2, 3, 4, 5, 6, 7, 8, 9]; G = 0; b = 0; yy = 0.0
for s in seeds:
    z = np.load(os.path.join(ROOT, "results", f"gram_synthetic_s{s}.npz")); G = G + z["G"]; b = b + z["b"]; yy += float(z["yy"])
keep = np.array(sorted(best["units"])); K = len(keep); assert K == best["K"]
g = np.linalg.solve(G[np.ix_(keep, keep)], b[keep]); r2 = 1 - (yy - b[keep] @ g) / yy
assert abs((1 - r2) - best["e2"]) < 1e-9, (1 - r2, best["e2"])
rows = pair[keep].ravel(); offs = sorted({(int(T["di"][r]), int(T["dj"][r])) for r in rows}); L = len(offs); assert L == best["loci"]
Ej = jonswap(); Yj = wrt(Ej).ravel(); Pj = g @ unit_responses(Ej)[0][keep]; c_chk = float((Pj * Yj).sum() / (Pj * Pj).sum()); shutil.rmtree(XNL)
name = f"swrt_{grid.TAG}_K{K}L{L}"; npz = os.path.join(HERE, f"{name}.npz"); w3 = os.path.join(ROOT, "ww3", "nl6", f"{name}.ww3op")
assert not os.path.exists(npz) and not os.path.exists(w3), "operator files exist already; not overwritten"
note = (f"SWRT (ecWAM grid {grid.TAG}: {grid.NA} directions x {grid.NK} frequencies from {grid.F0:.6g} Hz x 1.1, ecWAM cost basis of swrt_op36, 2026-11-04): the smallest ebar_2 on the 4x DIA line of cost_ec/time_sets_ec.json, "
        f"one_L({best['L']}, default_rng({best['seed']}), K = {K}) on {L} offsets; least-squares gains on the summed synthetic Gram, c = 1")
np.savez(npz, K=K, units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=c_chk, probes="synthetic", seeds=np.array(seeds),
         probe_fit_r2=r2, **{k: T[k][rows] for k in swrt.GEOM + ("amem", "node")}, **{k: T[k] for k in swrt.CONST}, note=note)
arr = eo.kernel_arrays(npz); o = np.lexsort((arr["dj"], arr["di"]))
for n in eo.NAMES: arr[n] = arr[n][o] if arr[n].ndim == 1 else np.asfortranarray(arr[n][:, o])
eo.write(w3, arr)
byKs = {}
for r in R: byKs.setdefault(r["K"], {})[r["seed"]] = dict(L=r["L"], e2=r["e2"], uses=r["loci"], cost=r["cost"], src=r.get("src"), seconds=r.get("seconds"))
cfg = {}
for u, gu in zip(keep, g):
    di, dj = int(T["di"][pair[u][0]]), abs(int(T["dj"][pair[u][0]])); cfg.setdefault(f"({di}, {dj})", []).append(round(float(gu), 3))
out = dict(note=note, rank=RANK, K=K, L=L, L_limit=best["L"], seed=best["seed"], e2=float(1 - r2), train_pct=float(100 * np.sqrt(1 - r2)), cost_dia=best["cost"],
           offsets=offs, units=[int(u) for u in keep], gains=[float(x) for x in g], gains_by_config=cfg, c_check=c_chk,
           npz=npz, ww3op=w3, rows=int(len(rows)), line=byKs)
json.dump(out, open(os.path.join(HERE, f"swrt_ec4x_K{KF}{f'_r{RANK}' if RANK else ''}.json" if KF else "swrt_ec4x.json"), "w"), indent=1)
print(f"{name}: K {K}, offsets {L} (limit {best['L']}, seed {best['seed']}), e2 {1 - r2:.5f} (training error {100 * np.sqrt(1 - r2):.2f} %), "
      f"cost {best['cost']:.3f}x DIA, gains {g.min():.3g}..{g.max():.3g} (negative {(g < 0).sum()}), c check {c_chk:.4f}")
print("offsets:", offs); print("wrote", npz, "and", w3)

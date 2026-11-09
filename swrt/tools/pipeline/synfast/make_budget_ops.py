"""Operator files of the NBEST distinct sets with the smallest e2 of a budget search (budget_<c>.json of search_budget.py), by the SWRT
recipe as make_k7_ops.py: least-squares gains on the summed synthetic Gram, c = 1, lines_ec/low_ops/<name>.npz and ww3/nl6/<name>.ww3op
(new files only); a set equal to an existing file of low_ops keeps that file's name, the others are swrt_ec_b<c>_K<K>L<L> on the
36 x 36 grid and swrt_<tag>_b<c>_K<K>L<L> on a grid of swrt_grids (TAG of its grid.py), with a letter added when another set of the same
K and L already has the name. Given set keys of budget_<c>.json instead of NBEST (as K8Llim7), it writes the files of those sets only.
usage: python3 make_budget_ops.py <grid root> <budget> [NBEST | <set key> ...]  ->  names on stdout (one line, in e2 order or in the
       order of the keys), synfast/budget_<c>_ops.json, or synfast/budget_<c>_keys.json for given keys (accumulated)"""
import os, sys, json, glob, numpy as np
ROOT, BUD = os.path.abspath(sys.argv[1]), float(sys.argv[2]); KEYS = [a for a in sys.argv[3:] if not a.isdigit()]
NB = int(sys.argv[3]) if len(sys.argv) > 3 and not KEYS else 12
sys.path[:0] = [ROOT, os.path.join(ROOT, "ww3", "nl6")]; os.chdir(ROOT)
import swrt, grid, export_operator as eo
from build_swrt import pair, T
HERE = os.path.join(ROOT, "synfast"); LOW = os.path.join(ROOT, "lines_ec", "low_ops"); NL6 = os.path.join(ROOT, "ww3", "nl6")
B = json.load(open(os.path.join(HERE, f"budget_{BUD:g}.json"))); tag = f"{BUD:g}".replace(".", "p"); os.makedirs(LOW, exist_ok=True)
GT = getattr(grid, "TAG", None); PRE = f"swrt_{GT}" if GT else "swrt_ec"; GRID = f"{grid.NA} x {grid.NK} x {grid.XFR}" + (f", {GT}" if GT else "")
OLD = {tuple(sorted(int(u) for u in np.load(f)["units"])): os.path.basename(f)[:-4] for f in glob.glob(os.path.join(LOW, "*.npz"))}
seen, names, res = set(), [], {}
for k, r in ([(k, B["sets"][k]) for k in KEYS] if KEYS else sorted(B["sets"].items(), key=lambda kv: kv[1]["e2"])):
    key = tuple(sorted(r["units"]))
    if key in seen: continue
    seen.add(key); K, L = r["K"], r["L"]; keep = np.array(key); g = np.array(r["g"])
    name = OLD.get(key)
    if name is None:                                                     # a new set: a letter is added while the name is taken
        name = f"{PRE}_b{tag}_K{K}L{L}"; i = 0
        while os.path.exists(os.path.join(LOW, name + ".npz")): name = f"{PRE}_b{tag}_K{K}L{L}{'abcdefghij'[i]}"; i += 1
    npz = os.path.join(LOW, name + ".npz"); w3 = os.path.join(NL6, name + ".ww3op")
    if not os.path.exists(npz):
        assert not os.path.exists(w3); rows = pair[keep].ravel()
        note = (f"SWRT for the speed budget {BUD:g} x DIA (ecWAM {GRID}, 2026-11-08): K = {K} on {L} offsets from search_budget.py "
                f"(e2 on the summed synthetic Gram, positive gains, limit L {r['L_limit']}, seed {r['seed']}), c = 1")
        np.savez(npz, K=K, units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=np.nan, probes="synthetic",
                 seeds=np.array([11, 1, 2, 3, 4, 5, 6, 7, 8, 9]), probe_fit_r2=1 - r["e2"],
                 **{q: T[q][rows] for q in swrt.GEOM + ("amem", "node")}, **{q: T[q] for q in swrt.CONST}, note=note)
        arr = eo.kernel_arrays(npz); o = np.lexsort((arr["dj"], arr["di"]))
        for n in eo.NAMES: arr[n] = arr[n][o] if arr[n].ndim == 1 else np.asfortranarray(arr[n][:, o])
        eo.write(w3, arr)
    res[name] = dict(K=K, L=L, cost=r["cost"], e2=r["e2"], units=list(key), g=r["g"], npz=npz, ww3op=w3,
                     offsets=sorted({(int(T["di"][q]), int(abs(T["dj"][q]))) for u in keep for q in pair[u]}))
    names.append(name)
    if len(names) == NB and not KEYS: break
OUT = os.path.join(HERE, f"budget_{BUD:g}_{'keys' if KEYS else 'ops'}.json")
if KEYS and os.path.exists(OUT): res = dict(json.load(open(OUT)), **res)
json.dump(res, open(OUT, "w"), indent=1)
for n in names: print(f"{n:24s} cost {res[n]['cost']:.2f}x e2 {res[n]['e2']:.4f} {res[n]['offsets']}", file=sys.stderr)
print(" ".join(names))

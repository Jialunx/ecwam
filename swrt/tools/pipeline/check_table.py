"""Check that the folded table, summed with WRT's run-time rules (swrt.Rows), is WRT; and what each rule matters.

  exact          all 1184 rows (both senses of di = 0), density filter on       -> should equal WRT
  ring_once      di = 0 rows of one sense only (each equal-frequency pair once)
  no_density     density filter off (qf_frac = 0)
  wide_window    add every pair the file stores outside WRT's run-time window (di > 4 or |dj| > 5); WRT never uses them
Spectra: the test set of testsets.py (single-point WAVEWATCH III scenario runs, by scenario family); reference = fresh WRT.
Also: the exact table on the test spectra whose peak lies less than five bins above the lowest frequency, where the grid no
longer resolves the spectrum (peak_below_bin5; none in the six-scenario test set), and on all the others (peak_ge_bin5).
usage: python3 check_table.py  ->  results/check_table.json"""
import os, json, numpy as np
import swrt
from bqf_decode import decode
from build_table import fold
from wrt_ref import consts
from testsets import load, GROUPS

HERE = os.path.dirname(os.path.abspath(__file__))
T = swrt.load_table(); E, Y, grp, _ = load()


def run(Rw, name, T0=None):
    P = np.concatenate([Rw.snl(E[s:s + 25]) for s in range(0, len(E), 25)])
    sk = 1 - ((P - Y) ** 2).sum((1, 2)) / (Y ** 2).sum((1, 2)); mr = np.abs(P - Y).max((1, 2)) / np.abs(Y).max((1, 2))
    out = {g: dict(n=int((grp == g).sum()), skill_median=float(np.median(sk[grp == g])), skill_p5=float(np.percentile(sk[grp == g], 5)),
                   skill_min=float(sk[grp == g].min()), maxrel_median=float(np.median(mr[grp == g])), maxrel_max=float(mr[grp == g].max()))
           for g in GROUPS}
    print(name, "  ".join(f"{g}: med {v['skill_median']:.10f} p5 {v['skill_p5']:.6f} min {v['skill_min']:.6f} maxrel {v['maxrel_median']:.1e}/{v['maxrel_max']:.1e}"
                          for g, v in out.items()), flush=True)
    return out, P


res = {}
res["exact"], P0 = run(swrt.Rows(T), "exact")
res["ring_once"], _ = run(swrt.Rows(T, np.flatnonzero(~((T["di"] == 0) & (T["dj"] < 0)))), "ring_once")
Tn = dict(T); Tn["qf_frac"] = np.array(0.0); Rn = swrt.Rows(Tn)
res["no_density"], Pn = run(Rn, "no_density")
rel = np.abs(Pn - P0) / np.abs(P0).max((1, 2), keepdims=True); sk0 = 1 - ((P0 - Y) ** 2).sum((1, 2)) / (Y ** 2).sum((1, 2))
skn = 1 - ((Pn - Y) ** 2).sum((1, 2)) / (Y ** 2).sum((1, 2))
res["no_density"]["all_cells_changed_over_1e-6"] = int((rel > 1e-6).sum()); res["no_density"]["all_cells_total"] = int(rel.size)
res["no_density"]["all_max_change_rel"] = float(rel.max()); res["no_density"]["all_max_skill_change"] = float(np.abs(skn - sk0).max())
print("density filter:", {k: v for k, v in res["no_density"].items() if k.startswith("all_")})
Q = decode(os.path.join(HERE, "quad99999.bqf")); C = consts()
from grid import NK, NKMEM, NAMEM                         # 36 x 50 copy
mk = [(m, a) for m in range(NKMEM) for a in range(NAMEM) if Q["nloc"][m, a] > 0]
jk = {m: (j1, j3) for j1 in range(1, NK + 1) for j3 in range(j1, NK + 1) for m in [j3 + (j1 - 1) * NK - j1 * (j1 - 1) // 2 - 1]}
dimax = max(jk[m][1] - jk[m][0] for m, a in mk); amax = max(a for m, a in mk)
Tw, _ = fold(Q, C, di_max=dimax, dj_max=amax); res["file_extent"] = dict(di_max=int(dimax), amem_max=int(amax + 1), wide_rows=len(Tw["di"]))
res["wide_window"], _ = run(swrt.Rows(Tw), f"wide_window (di <= {dimax}, amem <= {amax + 1}, {len(Tw['di'])} rows)")
low = np.argmax(E.sum(-1), 1) < 5                                              # peak less than five bins above the lowest frequency
_, scen = load()[2:]; full = swrt.Rows(T)
for key, m in (("peak_below_bin5", low), ("peak_ge_bin5", ~low)):
    if not m.any(): res[key] = dict(n=0, scenarios=[]); print(key, res[key], flush=True); continue   # no such spectrum in the test set
    Pm = np.concatenate([full.snl(E[m][s:s + 25]) for s in range(0, int(m.sum()), 25)])
    rl = np.abs(Pm - Y[m]).max((1, 2)) / np.abs(Y[m]).max((1, 2))
    res[key] = dict(n=int(m.sum()), scenarios=sorted(set(scen[m])), maxrel_median=float(np.median(rl)), maxrel_max=float(rl.max()))
    print(key, res[key], flush=True)
json.dump(res, open(os.path.join(HERE, "results", "check_table.json"), "w"), indent=1)

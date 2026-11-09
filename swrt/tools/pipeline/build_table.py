"""Fold WRT's working set in quad99999.bqf into the row table that sWRT selects from.

WRT's run-time window: k3/k1 <= 2.5 (di = jk3 - jk1 = 0..4) and |theta3 - theta1| <= 75 deg (amem - 1 = 0..5).
Deep-water loci on the log grid are self-similar, so one configuration (di, amem) has the same node offsets at
every jk1. Row = (config, node): geometry (u2, v2, u4, v4) from the reference bin jk1 = 20, weights zz(jk1) for
all 40 bins, matched by node index when the node count equals the reference count, otherwise by geometry
(tolerance 0.25 bin). A node whose k2 or k4 lies just beyond bin 40 at jk1 = 20 is stored there in WRT's out-of-grid
encoding (ik = 40, fraction (k - k40)/50), which is not its true offset; its geometry is taken from the nearest lower jk1
where it is inside the grid (at jk1 = 20 itself both corners clamp to bin 40 either way).
Direction senses exactly as WRT visits them (q_xnl4v4 loops ik3 = ik1..40, ia3 = 1..24):
  amem = 2..12: +dj and -dj (the -dj row is the mirror image, v -> -v);  amem = 1 (di > 0) and amem = 13: once;
  di = 0 is included with both senses, i.e. every equal-frequency pair is visited from both ends, as in WRT.
Also stored: WRT's run-time constants, so the sWRT kernel never needs the WRT library.
usage: python3 build_table.py  ->  data/table.npz, data/table.csv"""
import os, json, numpy as np
from bqf_decode import decode, offsets, raw_nodes, NK, NA
from wrt_ref import consts

from grid import JREF, DI_MAX, DJ_MAX                 # ecWAM copy: grid.py (reference bin 20, window di <= 4, |dj| <= 7)
HERE = os.path.dirname(os.path.abspath(__file__)); TOL = 0.25


def fold(Q, C, di_max=DI_MAX, dj_max=DJ_MAX):
    rows, st = [], dict(fallback_bins=0, unmatched_nodes=0, unmatched_zz=0.0, total_zz=0.0, edge_fixed=[])
    for di in range(0, di_max + 1):
        for dj in range(0, dj_max + 1):
            if di == 0 and dj == 0: continue                       # k1 = k3 is skipped by WRT
            jr = JREF if JREF + di <= NK else NK - di
            ref = offsets(Q, jr, jr + di, dj + 1); rr = raw_nodes(Q, jr, jr + di, dj + 1)
            for m in np.flatnonzero((rr[:, 0] >= NK) | (rr[:, 2] >= NK)) if len(rr) else []:
                for j in range(jr - 1, 0, -1):                       # nearest lower bin with the node inside the grid
                    r2 = raw_nodes(Q, j, j + di, dj + 1)
                    if len(r2) == len(rr) and r2[m, 0] < NK and r2[m, 2] < NK:
                        st["edge_fixed"].append((di, dj, int(m) + 1, float(ref[m, 0]), float(offsets(Q, j, j + di, dj + 1)[m, 0]), j))
                        ref[m, :4] = offsets(Q, j, j + di, dj + 1)[m, :4]; break
            W = np.zeros((len(ref), NK))
            for j1 in range(1, NK + 1 - di):
                nd = offsets(Q, j1, j1 + di, dj + 1); st["total_zz"] += nd[:, 4].sum() if len(nd) else 0.0
                if len(nd) == len(ref): W[:, j1 - 1] = nd[:, 4]; continue
                st["fallback_bins"] += 1
                for q in nd:                                        # match each stored node to its reference node
                    d = np.abs(ref[:, :4] - q[:4]).max(1) if len(ref) else np.array([9.0])
                    if d.min() < TOL: W[d.argmin(), j1 - 1] = q[4]
                    else: st["unmatched_nodes"] += 1; st["unmatched_zz"] += q[4]
            for m in range(len(ref)):
                for s in ((1,) if dj in (0, NA // 2) else (1, -1)):
                    rows.append((di, s * dj, dj + 1, m, ref[m, 0], s * ref[m, 1], ref[m, 2], s * ref[m, 3], W[m]))
    t = {k: np.array([r[i] for r in rows], int) for i, k in enumerate(("di", "dj", "amem", "node"))}
    t.update({k: np.array([r[i] for r in rows], float) for i, k in enumerate(("u2", "v2", "u4", "v4"), 4)})
    t["zz"] = np.array([r[8] for r in rows])
    t.update(freq=C["freq"], sigma=C["sigma"], k=C["k"], cg=C["cg"], dk=C["dk"], dtheta=C["dtheta"],
             a=C["k"] * C["dtheta"] * C["dk"],                        # area factor k dtheta dk
             b=C["k"] / C["cg"] * 2 * np.pi * C["sigma"],             # N(kx,ky) -> E(f,theta) factor
             scale=2 * C["q_scale"], qf_frac=C["qf_frac"], kpow=C["kpow"], depth_ref=C["depth"])
    return t, st


def mirror_units(t):
    """mirror partner (theta -> -theta about k1): for dj != 0 the same node with the opposite sense; for dj = 0 the node
    of the same configuration at (u2, -v2, u4, -v4) (dj = 0 loci are symmetric, nodes match to 1e-5 bin).
    A unit = {row, mirror}; sWRT selects units, so the operator keeps WRT's left/right symmetry."""
    di, dj, node = t["di"], t["dj"], t["node"]; n = len(di)
    key = {(a, b, c): i for i, (a, b, c) in enumerate(zip(di, dj, node))}; mirror = np.full(n, -1)
    for i in range(n):
        if dj[i] != 0: mirror[i] = key[(di[i], -dj[i], node[i])]; continue
        same = np.flatnonzero((di == di[i]) & (dj == 0)); G = np.column_stack([t[k][same] for k in ("u2", "v2", "u4", "v4")])
        d = np.abs(G - np.array([t["u2"][i], -t["v2"][i], t["u4"][i], -t["v4"][i]])).max(1); mirror[i] = same[d.argmin()]
        assert d.min() < 1e-4
    assert (mirror[mirror] == np.arange(n)).all() and (mirror != np.arange(n)).all()
    return mirror, np.unique(np.minimum(np.arange(n), mirror), return_inverse=True)[1]


if __name__ == "__main__":
    Q = decode(os.path.join(HERE, "quad99999.bqf")); C = consts()
    assert np.allclose(Q["sigma"], C["sigma"], rtol=1e-5)          # 36 x 50: Xnl's single-precision sigma is 2e-6 off 2 pi f at the top
    t, st = fold(Q, C); t["mirror"], t["unit"] = mirror_units(t)
    np.savez(os.path.join(HERE, "data", "table.npz"), **t,
             note="rows of WRT's own quadrature (quad99999.bqf), working window di 0..DI_MAX, |dj| <= DJ_MAX of grid.py; zz raw (N(kx,ky) units)")
    cols = ["di", "dj", "amem", "node", "u2", "v2", "u4", "v4"] + [f"zz_f{f:.4f}" for f in C["freq"]]
    with open(os.path.join(HERE, "data", "table.csv"), "w") as fh:
        fh.write("row," + ",".join(cols) + "\n")
        for i in range(len(t["di"])):
            fh.write(f"{i}," + ",".join(str(t[k][i]) for k in ("di", "dj", "amem", "node")) + "," + ",".join(f"{t[k][i]:.6f}" for k in ("u2", "v2", "u4", "v4"))
                     + "," + ",".join(f"{x:.6e}" for x in t["zz"][i]) + "\n")
    print(f"table: {len(t['di'])} rows, {len(set(zip(t['di'], t['dj'])))} (di, dj) configurations, {(t['dj'] == 0).sum()} rows with dj = 0, "
          f"mirror-pair units {t['unit'].max() + 1}")
    print(f"stored nodes in window: total zz {st['total_zz']:.6e}; bins matched by geometry {st['fallback_bins']}, "
          f"unmatched nodes {st['unmatched_nodes']} (zz share {st['unmatched_zz'] / st['total_zz']:.2e})")
    print("edge-encoded reference nodes re-read (di, dj, node, u2 at jk1=20, u2 used, from jk1):", st["edge_fixed"])
    dm = [float(max(abs(t["u2"][i] - t["u2"][t["mirror"][i]]), abs(t["v2"][i] + t["v2"][t["mirror"][i]]),
                    abs(t["u4"][i] - t["u4"][t["mirror"][i]]), abs(t["v4"][i] + t["v4"][t["mirror"][i]]))) for i in np.flatnonzero(t["dj"] == 0)]
    json.dump(dict(rows=len(t["di"]), configs=len(set(zip(t["di"], t["dj"]))), rows_dj0=int((t["dj"] == 0).sum()), units=int(t["unit"].max() + 1),
                   fallback_bins=st["fallback_bins"], unmatched_nodes=st["unmatched_nodes"], unmatched_zz_share=st["unmatched_zz"] / st["total_zz"],
                   edge_fixed=st["edge_fixed"], dj0_mirror_position_mismatch_max=max(dm)),
              open(os.path.join(HERE, "results", "table_build.json"), "w"), indent=1)

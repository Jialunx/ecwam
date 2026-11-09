"""Write an sWRT operator (operators/swrt_K*.npz, from build_swrt.py) as the plain-text file read by WAVEWATCH III's NL6
(ww3/WW3/model/src/w3snl6md.F90). The file holds the row arrays of the compiled kernel (bench/kernels.f90, swrt_batch),
built as bench/bench.py rows_of(d, g_row) builds them, so that W3SNL6 only reads and applies them.

usage: python3 export_operator.py                          ->  swrt_K30.ww3op (next to this script, from operators/swrt_K30.npz)
       python3 export_operator.py --op path/to/swrt_K50.npz  ->  swrt_K50.ww3op
       python3 export_operator.py --out /some/run/dir/swrt_operator.ww3op

Format (version 1, plain text, read with list-directed Fortran READ; integers %d, reals %.9e, which round-trips REAL*4):
    WW3-SWRT-OPERATOR          magic line
    1                          format version
    nk na nr K                 frequencies, directions, rows, mirror-pair units (nr = 2 K for build_swrt.py operators)
    swrt_K30.npz               source operator file (name only)
followed by one block per array: a line with the array name, then its values (fastest index first):
    sigma(nk)                  radian frequencies of the table grid (rad/s); W3SNL6 checks WAVEWATCH III's SIG against them
    di dj(nr)                  k3 offset in frequency / direction bins
    i2 f2 j2 g2(nr)            floor and fraction of the k2 offset: frequency (i2, f2) and direction (j2, g2); i4 f4 j4 g4 for k4
    w fa fp(nk, nr)            zz * gain * scale (0 where k3 is off the grid); a(i3) b(i1); a(i1) b(i3)
    conv(nk)                   E(f, theta) -> N(kx, ky): cg / (2 pi sigma k)
    qfrac                      density filter: skip a pair if N1 and N3 < qfrac * max(N) * kpow(i3)
    kpow(nk)                   (k1 / k)^7.5
The gains are the operator's g_row times its c (c = 1 for every operator written by build_swrt.py, as used by swrt.SWRT).
The rows are written grouped by resonance curve (sorted by di, then dj; stable), because W3SNL6's default kernel shares
the k1 and k3 work of the rows of one curve. The order of the rows changes the result only by round-off."""
import os, sys, argparse, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); sys.path.insert(0, ROOT)
import swrt

MAGIC, VERSION, TPI = "WW3-SWRT-OPERATOR", 1, 2 * np.pi
NAMES = ("di", "dj", "i2", "f2", "j2", "g2", "i4", "f4", "j4", "g4", "w", "fa", "fp")


def rows_of(d, g_row, T):
    """bench/bench.py rows_of, unchanged except that the run-time constants T are an argument (bench uses data/table.npz)"""
    NK = len(T["sigma"])
    R = len(d["di"]); di = np.asarray(d["di"], int); I = np.arange(NK)
    i3 = I[None] + di[:, None]; ok = i3 <= NK - 1; i3c = np.minimum(i3, NK - 1)
    w = np.where(ok, np.asarray(d["zz"], float) * np.asarray(g_row, float)[:, None] * float(T["scale"]), 0.0)
    fa, fp = T["a"][i3c] * T["b"][None, :], T["a"][None, :] * T["b"][i3c]
    fl = lambda x: np.floor(np.asarray(x, float)).astype(np.int32); fr = lambda x: (np.asarray(x, float) - np.floor(np.asarray(x, float))).astype(np.float32)
    F = lambda x: np.asfortranarray(x.T, dtype=np.float32)                     # (nk, R)
    return (di.astype(np.int32), np.asarray(d["dj"], np.int32), fl(d["u2"]), fr(d["u2"]), fl(d["v2"]), fr(d["v2"]),
            fl(d["u4"]), fr(d["u4"]), fl(d["v4"]), fr(d["v4"]), F(w), F(fa), F(fp))


def kernel_arrays(op):
    """all swrt_batch inputs of an operator file, as bench.py builds them, from the operator's own constants"""
    d = np.load(op); T = {k: np.asarray(d[k]) for k in swrt.CONST}
    Tt = swrt.load_table()
    diff = [k for k in swrt.CONST if not np.array_equal(T[k], np.asarray(Tt[k]))]
    if diff: print(f"note: {os.path.basename(op)} constants differ from data/table.npz in {diff}; using the operator's own")
    sig, k, cg = T["sigma"], T["k"], T["cg"]
    out = dict(zip(NAMES, rows_of(d, np.asarray(d["g_row"], float) * float(d["c"]), T)))
    out.update(sigma=np.asarray(sig, float), conv=(cg / (TPI * sig * k)).astype(np.float32), qfrac=np.float32(T["qf_frac"]),
               kpow=T["kpow"].astype(np.float32), K=int(d["K"]), nk=len(sig), na=swrt.NA, nr=len(d["di"]), source=os.path.basename(op))
    return out


def write(path, a, per_line=8):
    def block(f, name, v, fmt):
        f.write(name + "\n"); v = np.ravel(v, order="F")
        for s in range(0, len(v), per_line): f.write(" ".join(fmt % x for x in v[s:s + per_line]) + "\n")
    with open(path, "w") as f:
        f.write(f"{MAGIC}\n{VERSION}\n{a['nk']} {a['na']} {a['nr']} {a['K']}\n{a['source']}\n")
        block(f, "sigma", a["sigma"], "%.9e")
        for n in NAMES:
            v = a[n]; assert np.all(np.isfinite(v)), n
            block(f, n, v, "%d" if v.dtype.kind == "i" else "%.9e")
        block(f, "conv", a["conv"], "%.9e"); block(f, "qfrac", [a["qfrac"]], "%.9e"); block(f, "kpow", a["kpow"], "%.9e")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--op", default=os.path.join(ROOT, "operators", "swrt_K30.npz"), help="operator .npz (default operators/swrt_K30.npz)")
    ap.add_argument("--out", default=None, help="output file (default <op name>.ww3op next to this script)")
    a = ap.parse_args()
    out = a.out or os.path.join(HERE, os.path.splitext(os.path.basename(a.op))[0] + ".ww3op")
    arr = kernel_arrays(a.op); order = np.lexsort((arr["dj"], arr["di"]))                 # rows of a curve together
    for n in NAMES: arr[n] = arr[n][order] if arr[n].ndim == 1 else np.asfortranarray(arr[n][:, order])
    write(out, arr)
    print(f"wrote {out}: K = {arr['K']}, {arr['nr']} rows, nk = {arr['nk']}, na = {arr['na']} (from {arr['source']})")

"""Integration filter of SWRT files on training spectra only (from swrt_op36/synfast/dyn_syn.py, 2026-11-08). A set whose transfer
blows up in time is dropped: from training spectra, dE/dt = S_nl(E) is integrated with the file, and the file passes when the run stays
finite. Tools copy: only this finite / blows-up test is kept. The trajectory score J_dyn, and with it the WRT trajectories of
dyn_ref.npz, are left out; the choice among the finite sets is the smallest training e2 (make_swrt_files.sh).
  spectra   NSP per group of practice draw SPEC_SEED (in the training pool 11, 1-9; narrow single-peak = directional spread below 20 deg,
            broad single-peak, two-peak; drawn with seed 7). The held-out draw 101 and the model tests are not used.
  steps     semi-implicit as WAM and WAVEWATCH III, dE = dt S / (1 - dt min(D, 0)) with D = dS/dE on the diagonal (Rows.diag), E >= 0,
            NSTEP steps of dt_p = CFL / r_p per spectrum, r_p = sum |S_WRT| w / sum E w at t = 0 (WRT = the full table, = WRT to 1e-5),
            the same dt as the original's WRT run
  finite    every checkpoint (every EVERY steps) of every spectrum finite
usage: python3 dyn_syn.py <grid root> <npz file> [...]  ->  <grid root>/synfast/dyn_syn.json (accumulated by file name)"""
import os, sys, json, time, fcntl, numpy as np
from multiprocessing import Pool
ROOT = os.path.abspath(sys.argv[1]); FILES = [os.path.abspath(f) for f in sys.argv[2:]]; sys.path.insert(0, ROOT); os.chdir(ROOT)
import build_swrt as bs, swrt
OD = os.path.join(ROOT, "synfast"); OUT = os.path.join(OD, "dyn_syn.json")
SPEC_SEED, NSP, NSTEP, CFL, EVERY = 1, 12, 100, 0.05, 10
T = bs.T; NA = int(bs.DIRS.size); dth = 2 * np.pi / NA
w = (T["cg"] * (T["a"] / (T["k"] * dth)) / (2 * np.pi) * dth)[:, None] * np.ones((1, NA)); th = np.deg2rad(bs.DIRS)
Z = np.load(os.path.join(ROOT, "data", f"synthetic_s{SPEC_SEED}.npz")); EA, lab = Z["E"].astype(float), Z["lab"]
Ed = (EA * w[None]).sum(1); spr = np.degrees(np.sqrt(2 * (1 - np.hypot((Ed * np.cos(th)).sum(1), (Ed * np.sin(th)).sum(1)) / Ed.sum(1))))
rng = np.random.default_rng(7); GRP = {}
for g, m in (("narrow", (lab == 0) & (spr < 20)), ("broad", (lab == 0) & (spr >= 20)), ("two", lab == 1)):
    GRP[g] = np.sort(rng.choice(np.flatnonzero(m), NSP, replace=False))
IDX = np.concatenate(list(GRP.values())); E0 = EA[IDX]


def run(snl, dia, E, dt):
    out = [E.copy()]
    for n in range(1, NSTEP + 1):
        S, D = snl(E), dia(E)
        E = np.maximum(E + dt[:, None, None] * S / (1 - dt[:, None, None] * np.minimum(D, 0.0)), 0.0)
        if n % EVERY == 0: out.append(E.copy())
    return np.stack(out, 1)                                             # (spectrum, checkpoint, frequency, direction)


def dt_chunk(i):                                                        # the original's ref_chunk without the WRT run
    E = E0[i:i + 3]; S = bs.R.snl(E); return CFL * (E * w[None]).sum((1, 2)) / (np.abs(S) * w[None]).sum((1, 2))


def cand(a):
    f, dt = a; op = swrt.SWRT(f); t0 = time.time()
    XC = run(op.snl, op.diag, E0.copy(), dt)
    return os.path.basename(f)[:-4], dict(finite=bool(np.isfinite(XC).all()), file=f, seconds=time.time() - t0)


if __name__ == "__main__":
    dt = np.concatenate([dt_chunk(i) for i in range(0, len(E0), 3)])
    def save(name, r):                                                   # merge under a lock
        with open(OUT + ".lock", "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX); res = json.load(open(OUT)) if os.path.exists(OUT) else {}; res[name] = r
            json.dump(res, open(OUT + ".tmp", "w"), indent=1); os.replace(OUT + ".tmp", OUT)
    if FILES:
        with Pool(min(int(os.environ.get("DYN_WORKERS", 4)), len(FILES))) as p:
            for name, r in p.imap_unordered(cand, [(f, dt) for f in FILES]):
                save(name, r)
                print(f"{name:24s} {'finite' if r['finite'] else 'blows up'} ({r['seconds']:.0f} s)", flush=True)

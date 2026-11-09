"""Synthetic practice (probe) spectra for step 9: a JONSWAP family generated from formulas and one random seed; no model data.
  900 single-peak spectra: fp log-uniform 0.05-0.30 Hz, gamma uniform 1-7, directional spread at the peak log-uniform 4-60 deg
      (cos^2s(dtheta/2) with s solved on the 36-bin grid), random mean direction, Hs 2 m;
  600 two-peak spectra: swell fp1 log-uniform 0.05-0.12 Hz, gamma 1-7, spread 4-20 deg; wind sea fp2 = (1.3-3.0) fp1, gamma 1-5,
      spread 20-50 deg; independent random directions; swell scaled by 0.2-5 relative to the wind sea (each Hs 1 m before scaling).
Same family as nlrt/lrt_build_synthetic.py, without WRT labels, but with the direction difference wrapped: the nlrt version clipped
cos((theta - theta_p)/2) at zero without wrapping, which cut the distribution on one side when theta_p was near 360 deg (fixed 2026-09-23),
so seed 11 no longer reproduces nlrt/data/zeroshot_corpus_synthetic.npz.
usage: python3 make_synthetic.py [seed ...]  ->  data/synthetic_s{seed}.npz   (default seed 11)"""
import os, sys, numpy as np

from grid import NA as ND, FREQ, DIRS
HERE = os.path.dirname(os.path.abspath(__file__)); G0 = 9.81
f = FREQ.copy(); df = np.gradient(f); th = np.deg2rad(DIRS)   # 36 x 50 copy: bin centres (as the Grams of candidates/grid36x50)


def spread_deg(D):
    D = D / D.sum(); C, S = (D * np.cos(th)).sum(), (D * np.sin(th)).sum(); return np.degrees(np.sqrt(2 * (1 - np.hypot(C, S))))


def dirdist(theta_p, target):
    """cos^(2s)((theta - theta_p)/2) on the 36-bin grid, s solved so that the discrete spread equals target (deg)"""
    base = np.abs(np.cos(0.5 * (th - np.deg2rad(theta_p)))); lo, hi = 0.5, 4000.0      # |cos| = cos of the wrapped half-difference
    for _ in range(60):
        s = np.sqrt(lo * hi); D = base ** (2 * s)
        lo, hi = (s, hi) if spread_deg(D + 1e-300) > target else (lo, s)
    D = base ** (2 * s); return D / (D.sum() * 2 * np.pi / ND)


def jonswap1d(fp, gamma, Hs):
    sig = np.where(f <= fp, 0.07, 0.09); rr = np.exp(-((f - fp) ** 2) / (2 * (sig * fp) ** 2))
    S1 = G0 ** 2 / ((2 * np.pi) ** 4 * f ** 5) * np.exp(-1.25 * (fp / f) ** 4) * gamma ** rr
    return S1 * (Hs / 4) ** 2 / (S1 * df).sum()


def make(seed):
    rng = np.random.default_rng(seed); E, lab = [], []
    for _ in range(900):                                   # single peak
        fp = np.exp(rng.uniform(np.log(0.05), np.log(0.30))); gam = rng.uniform(1.0, 7.0); sp = np.exp(rng.uniform(np.log(4.0), np.log(60.0)))
        E.append(jonswap1d(fp, gam, 2.0)[:, None] * dirdist(rng.uniform(0, 360), sp)[None, :]); lab.append(0)
    for _ in range(600):                                   # swell + wind sea
        fp1 = np.exp(rng.uniform(np.log(0.05), np.log(0.12))); fp2 = fp1 * rng.uniform(1.3, 3.0)
        e1 = jonswap1d(fp1, rng.uniform(1, 7), 1.0)[:, None] * dirdist(rng.uniform(0, 360), np.exp(rng.uniform(np.log(4.0), np.log(20.0))))[None, :]
        e2 = jonswap1d(fp2, rng.uniform(1, 5), 1.0)[:, None] * dirdist(rng.uniform(0, 360), rng.uniform(20.0, 50.0))[None, :]
        E.append(e1 * rng.uniform(0.2, 5.0) + e2); lab.append(1)
    return np.stack(E).astype(np.float32), np.array(lab)


if __name__ == "__main__":
    for seed in [int(a) for a in sys.argv[1:]] or [11]:
        E, lab = make(seed)
        np.savez(os.path.join(HERE, "data", f"synthetic_s{seed}.npz"), E=E, lab=lab, seed=seed, note=__doc__.split("\n")[0])
        print(f"seed {seed}: {len(E)} spectra ({(lab == 0).sum()} single-peak, {(lab == 1).sum()} two-peak)")

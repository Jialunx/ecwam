"""Test spectra of the ecWAM copy (replaces swrt_g36/testsets.py, which loads WAVEWATCH III single-point runs that do not exist on this
grid): 300 held-out synthetic spectra from make_synthetic.py seed 101 (not in the practice pool 11, 1-9), the first 150 single-peak
(rows 0-149) and the first 150 two-peak (rows 900-1049), each with a fresh WRT reference on this grid (as keydoc/swrt/candidates/
grid36x50 step 1 / step 5 did). Used only to check the table and to report the operator's error, never to choose anything.
    python3 testsets.py         ->  data/test_syn_s101.npz  (WRT on 300 spectra)
    from testsets import load, GROUPS
    E, Y, grp, scen = load()    # E, Y (n, NK, NA); groups 'single' / 'two-peak'; scen = the same labels"""
import os, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PATHS = {"test": os.path.join(HERE, "data", "test_syn_s101.npz")}
GROUPS = ("single", "two-peak")
SEED, ROWS = 101, np.r_[0:150, 900:1050]


def load(name="test"):
    z = np.load(PATHS[name])
    return z["E"].astype(float), z["Y"].astype(float), z["family"].astype(str), z["scenario"].astype(str)


if __name__ == "__main__":
    import time
    from wrt_ref import wrt
    z = np.load(os.path.join(HERE, "data", f"synthetic_s{SEED}.npz")); E = z["E"][ROWS].astype(float)
    fam = np.where(z["lab"][ROWS] == 0, "single", "two-peak"); t0 = time.time(); Y = wrt(E)
    np.savez(PATHS["test"], E=E, Y=Y, family=fam, scenario=fam, rows=ROWS, seed=SEED, note=__doc__.split("\n")[0])
    print(f"{len(E)} spectra, WRT {time.time() - t0:.0f} s -> {PATHS['test']}")

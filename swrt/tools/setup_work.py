"""Work folder of make_swrt_files.sh for one spectral grid of ecWAM (from swrt_grids/setup_grid.py, 2026-11-04): a copy of pipeline/
(the scripts of swrt_grids/<tag> and swrt_op36/synfast, adapted as their headers say), cost_ec/time_sets_ec.json (ecWAM's cost model
of the 36 x 36 grid, used for every grid) and grid.py. Code only; the data, Grams and searches are made by make_swrt_files.sh.
The frequency ratio of ecWAM is fixed (FRATIO = 1.1, yowfred.F90), so XFR = 1.1. Window of WRT on the grid (Xnl: k3/k1 <= 2.5,
|theta3 - theta1| <= 75 deg): di <= 4 (1.1^8 = 2.14 <= 2.5 < 1.1^10), |dj| <= floor(75 / dtheta). JREF is the 1-based bin nearest
0.21 Hz, as in swrt_op36. The tag is d<NANG>[_n<NFRE>]_f<first frequency>, as the tags of swrt_grids (d24_f0418, d24_n40_f0345).
An existing work folder is reused when its grid.py is the same grid (a resumed run).
usage: python3 setup_work.py <NANG> <NFRE> <fr1> <ifre1> <work folder>   ->  the tag on stdout"""
import os, sys, shutil, math

HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, "pipeline")


def tag_of(na, nk, ifre1, fr1):
    f0 = fr1 / 1.1 ** (ifre1 - 1); f = f"{f0:.4f}"[2:].rstrip("0")                 # 0.0418 -> 0418, 0.0550 -> 055
    return f"d{na}" + (f"_n{nk}" if nk != 36 else "") + f"_f{f}"


def grid_py(tag, na, nk, ifre1, fr1):
    f0 = fr1 / 1.1 ** (ifre1 - 1); dth = 360.0 / na; djm = int(math.floor(75.0 / dth + 1e-9))
    jref = min(range(1, nk + 1), key=lambda j: abs(f0 * 1.1 ** (j - 1) - 0.21))
    assert jref + 4 <= nk, f"JREF {jref} + 4 > NFRE {nk}: the grid ends too close to 0.21 Hz"
    return f'''"""The spectral grid of this work folder ({tag}, ecwam/swrt/tools/setup_work.py): {na} directions x {nk} frequencies.
NFRE = {nk} source-term frequencies F0 x 1.1^i, F0 = fr1 / 1.1^(ifre1 - 1) = {fr1!r} / 1.1^{ifre1 - 1} = {f0:.10g} Hz, i = 0..{nk - 1}
({f0:.4g}-{f0 * 1.1 ** (nk - 1):.3g} Hz), {na} directions ({dth:g} deg), deep water. Every module of this copy takes its grid from here.
WRT's window (Xnl qf_krat = 2.5, qf_dmax = 75 deg): di <= 4, |dj| <= {djm} ({djm * dth:g} deg). JREF: the bin nearest 0.21 Hz (1-based), {jref}
({f0 * 1.1 ** (jref - 1):.4f} Hz)."""
import numpy as np

NK, NA, XFR = {nk}, {na}, 1.1
IFRE1, FR1 = {ifre1}, {fr1!r}                        # ecWAM settings of this grid (FR(IFRE1) = FR1)
F0 = FR1 / XFR ** (IFRE1 - 1)                        # ecWAM FR(1)
DTH_DEG = 360.0 / NA
FREQ = F0 * XFR ** np.arange(NK)                    # Hz
DIRS = DTH_DEG / 2 + DTH_DEG * np.arange(NA)        # deg, bin centres (only the spacing matters to WRT)
DI_MAX, DJ_MAX = 4, {djm}                               # WRT's run-time window on this grid
JREF = {jref}                                           # anchor bin whose locus geometry defines the table rows
NKMEM, NAMEM = NK * (NK + 1) // 2, NA // 2 + 1      # frequency pairs (jk3 >= jk1) and direction differences stored by Xnl
TAG = "{tag}"
'''


def main(na, nk, fr1, ifre1, D):
    tag = tag_of(na, nk, ifre1, fr1); g = grid_py(tag, na, nk, ifre1, fr1); gp = os.path.join(D, "grid.py")
    if os.path.exists(gp):
        assert open(gp).read() == g, f"{D} holds another grid ({gp}); use a new work folder"
        return tag
    for d, _, fs in os.walk(SRC):
        for f in fs:
            if f.endswith((".py", ".sh")):
                o = os.path.join(D, os.path.relpath(os.path.join(d, f), SRC)); os.makedirs(os.path.dirname(o), exist_ok=True); shutil.copy2(os.path.join(d, f), o)
    os.makedirs(os.path.join(D, "cost_ec"), exist_ok=True); shutil.copy2(os.path.join(HERE, "time_sets_ec.json"), os.path.join(D, "cost_ec"))
    for d in ("data", "results", "operators", "ext", "files"): os.makedirs(os.path.join(D, d), exist_ok=True)
    open(gp, "w").write(g)
    return tag


if __name__ == "__main__":
    na, nk, fr1, ifre1, D = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), os.path.abspath(sys.argv[5])
    assert na > 0 and nk > 4 and ifre1 >= 1, "NANG > 0, NFRE > 4, ifre1 >= 1"
    os.makedirs(D, exist_ok=True); print(main(na, nk, fr1, ifre1, D))

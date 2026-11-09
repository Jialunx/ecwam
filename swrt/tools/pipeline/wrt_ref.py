"""The WRT reference (Van Vledder's GurboQuad / xnl4v5, f2py build in ext/xnl_wrt*.so) and the grid it runs on.

Configuration used for every label in this folder (ecWAM copy): 36 frequencies 0.0345227107 Hz (ecWAM FR(1)) x 1.1^i, 36 directions (10 deg),
iquad = 2 (deep-water loci, Webb coupling, WAM depth scaling), iq_grid = 3 (full circle), depth 1000 m,
f^-5 tail. WRT reads quad99999.bqf from the working directory, so every call runs inside this folder.

    from wrt_ref import wrt, consts
    S = wrt(E)          # E (40,24) or (n,40,24) energy density E(f,theta) [m^2/Hz/rad] -> Snl, same layout
"""
import os, sys, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "ext"))
from xnl_wrt import m_constants, m_xnldata as xn, serv_xnl4v5 as sv   # noqa: E402

from grid import NK, NA, FREQ, DIRS                    # ecWAM copy: 0.03453 Hz x 1.1^i, 10-deg bins
TPI, G, DEPTH = 2 * np.pi, 9.81, 1000.0
SIG = (TPI * FREQ).astype(np.float32)
ANG = (np.arange(NA) * TPI / NA).astype(np.float32)
_ready = False


@contextlib.contextmanager
def _here():
    cwd = os.getcwd(); os.chdir(HERE)
    try: yield
    finally: os.chdir(cwd)


def _init():
    global _ready
    if _ready: return
    with _here():
        m_constants.init_constants()
        xn.xnl_init(SIG, DIRS.astype(np.float32), np.float32(-5.0), np.float32(G), np.array([DEPTH], np.float32), 2, 3, 0)
    _ready = True


def wrt(E):
    """exact WRT transfer dE/dt(f, theta) for energy spectra E(f, theta)"""
    _init(); E = np.asarray(E, float); one = E.ndim == 2; E = E[None] if one else E
    out = np.empty_like(E)
    with _here():
        for n, e in enumerate(E):
            A = (e / (TPI * SIG[:, None])).astype(np.float32)          # A(sigma, theta) = E(f, theta) / (2 pi sigma)
            x, _, _ = xn.xnl_main(A, SIG, ANG, np.float32(DEPTH), 2, 0)
            out[n] = x * (TPI * SIG[:, None])
    return out[0] if one else out


def consts():
    """WRT's run-time constants, read from the module after one call (q_scale is only set inside q_xnl4v4)"""
    _init(); wrt(np.outer(np.exp(-((FREQ - 0.1) / 0.02) ** 2), np.cos(np.deg2rad(DIRS - 180) / 2) ** 2))
    k = np.array([sv.z_wnumb(s, np.float32(DEPTH), np.float32(G)) for s in SIG], float)
    cg = np.array([sv.z_cmpcg(s, np.float32(DEPTH), np.float32(G)) for s in SIG], float)
    dk = np.asarray(sv.z_steps(k.astype(np.float32)), float)
    return dict(freq=FREQ, sigma=SIG.astype(float), k=k, cg=cg, dk=dk, dtheta=float(xn.q_delta), q_scale=float(xn.q_scale),
                qf_frac=float(xn.qf_frac), kpow=(k[0] / k) ** 7.5, iq_dscale=int(xn.iq_dscale), depth=DEPTH, g=G)


def dscale(E, depth=DEPTH):
    """WAM depth-scaling factor that xnl_main applies after the integral (q_dscale), per spectrum"""
    E = np.asarray(E, float); E = E[None] if E.ndim == 2 else E
    s = SIG.astype(float); dsig = np.asarray(sv.z_steps(SIG), float); kk = np.array([sv.z_wnumb(x, np.float32(depth), np.float32(G)) for x in SIG], float)
    A = E / (TPI * s[None, :, None]); dn = A * dsig[None, :, None] * (TPI / NA)
    s0 = dn.sum((1, 2)); sk = (dn / np.sqrt(kk)[None, :, None]).sum((1, 2))
    kms = (s0 / np.where(sk > 0, sk, 1)) ** 2; kd = np.maximum(0.5, 0.75 * kms * depth)
    return np.where(s0 > 0, 1 + 5.5 / kd * (1 - 5 / 6 * kd) * np.exp(-1.25 * kd), 1.0)

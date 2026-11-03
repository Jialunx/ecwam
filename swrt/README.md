# SWRT for nonlinear four-wave interactions in ecWAM: installation and quick guide

SWRT (sparse Webb–Resio–Tracy) computes the nonlinear four-wave interactions S_nl from K configurations of the WRT
quadrature table, each a quadruplet and its mirror image, with one optimized coefficient per configuration. This patch adds
SWRT to ecWAM as an alternative to the Discrete Interaction Approximation (DIA) in `SNONLIN`. In a 24-h global test, SWRT
keeps the significant wave height within 1 % of WRT, compared with 15 % for DIA, and takes 1.6 times the time per step of
DIA.

The changes, the SWRT files and this guide are on the branch `swrt` of
[Jialunx/ecwam](https://github.com/Jialunx/ecwam/tree/swrt), a fork of [ecmwf-ifs/ecwam](https://github.com/ecmwf-ifs/ecwam)
at commit `e1779dc0` (ecWAM 1.5.14).

## 1. Changes to ecWAM

The patch adds one source file and changes eight others. The default `INLSCHEME = 0` leaves ecWAM unchanged.

| File | Change |
|---|---|
| `src/ecwam/snonlin_swrt.F90` (new) | `INISNONLIN_SWRT` reads the SWRT file and checks it against the spectral grid. `SNONLIN_SWRT` takes the arguments of `SNONLIN` and adds S_nl to `SL` and its diagonal term to `FLD`. |
| `src/ecwam/yowstat.F90`, `mpuserin.F90` | New switch `INLSCHEME` in namelist `NALINE`, 0 for DIA (default) and 1 for SWRT. |
| `src/ecwam/userin.F90` | Writes the scheme to `wam.log` and stops for other values. |
| `src/ecwam/initmdl.F90` | Calls `INISNONLIN_SWRT` when `INLSCHEME = 1`. |
| `src/ecwam/implsch.F90`, `wdfluxes.F90` | Call `SNONLIN_SWRT` instead of `SNONLIN` when `INLSCHEME = 1`. |
| `src/ecwam/CMakeLists.txt` | Adds the new file to the build. |
| `share/ecwam/scripts/ecwam_run_model.sh` | Passes `isnonlin`, `inlscheme` and `swrt_operator` of the yml settings to the namelist and the environment. |

Four properties of `SNONLIN_SWRT` matter for ecWAM.

- **Vectorization.** The innermost loop runs over the grid points of a block (NPROMA), for SWRT and WRT alike. The
  arithmetic is in single precision, and subnormal numbers are flushed to zero inside this routine only.
- **Directions.** The anticlockwise direction index of the SWRT file is mapped onto the clockwise index of ecWAM.
- **Depth.** The SWRT files hold the deep-water interactions, and the shallow-water factor of `SNONLIN` (`ISNONLIN`) scales
  them as it scales DIA.
- **Diagonal term.** ∂S_nl/∂F is computed as in Xnl and enters the implicit integration of ecWAM.

## 2. Installation

The branch can be cloned directly.

```bash
git clone -b swrt https://github.com/Jialunx/ecwam
```

An existing clone of ecWAM at commit `e1779dc0` takes the same changes from the patch.

```bash
git checkout e1779dc0aec1210ab6f5339a2582fbb7b6b18e6b
curl -LO https://raw.githubusercontent.com/Jialunx/ecwam/swrt/swrt/ecwam_swrt_e1779dc0.patch
git apply ecwam_swrt_e1779dc0.patch
```

ecWAM then builds as before and needs no new library. The SWRT files are in the folder `swrt/files` of the branch and go
to any directory the model can read.

```bash
curl -LO https://raw.githubusercontent.com/Jialunx/ecwam/swrt/swrt/files/swrt_ec4x_K31L13.ww3op
curl -LO https://raw.githubusercontent.com/Jialunx/ecwam/swrt/swrt/files/wrt_table_K807.ww3op
```

| SWRT file | Method | Configurations | Offsets | Use |
|---|---|---|---|---|
| `swrt_ec4x_K31L13.ww3op` | SWRT | 31 | 13 | runs with SWRT, selected for four times the S_nl cost of DIA |
| `wrt_table_K807.ww3op` | WRT | 807 | 74 | reference runs, all configurations with coefficient 1 |

## 3. Spectral grid

The SWRT files depend on the spectral grid only, not on the spatial grid or the time step. They are built for the spectral
grid of the ECMWF operational configuration, with 36 directions and 36 source-term frequencies from 0.0345 Hz with a ratio
of 1.1 (fr1 = 4.177248E-02 at ifre1 = 3). The same files therefore serve any spatial resolution with this spectral grid, such
as O48 or O1280.

A change in the number of directions, the number of frequencies, the first frequency or the frequency ratio needs a new SWRT
file. The model compares the file with its spectral grid at the start and stops with a message when they differ. The scripts
in `tools/` of the SWRT code archive ([doi:10.5281/zenodo.23654016](https://doi.org/10.5281/zenodo.23654016)) build a new
file from Xnl for any spectral grid in about an hour, and the code of ecWAM stays unchanged. The first frequency of the new
file is the first source-term frequency of ecWAM, fr1 / 1.1^(ifre1−1).

## 4. Running

| Setting | Where | Value |
|---|---|---|
| `INLSCHEME` | namelist `NALINE` | 1 for SWRT, 0 for DIA (default) |
| `WAM_SWRT_OPERATOR` | environment | path of the SWRT file |
| `ISNONLIN` | namelist `NALINE` | shallow-water scaling 0, 1 or 2, as for DIA |
| `WAM_SWRT_DIAG` | environment | `0` switches the diagonal term off (optional) |

The standalone scripts read the same settings from three lines of the yml file.

```yaml
physics:
    inlscheme: 1
    isnonlin: 0
    swrt_operator: /abs/path/swrt_ec4x_K31L13.ww3op
```

The log file `wam.log` names the scheme, the SWRT file, and its numbers of configurations and offsets. The model stops with a
message when the variable is unset or the file is missing.

## 5. Test results

### 5.1 Synthetic spectra

The error ε of the direction-integrated S_nl (van Vledder 2006, Eq. 58) is computed against WRT on 300 synthetic spectra
that were not used in the selection, half of them single-peaked and half bimodal. SWRT follows WRT closely, whereas the
error of DIA exceeds the magnitude of S_nl itself.

| Method | Median ε (%) | 95th percentile of ε (%) |
|---|---|---|
| SWRT with K = 31 | 8 | 18 |
| DIA of ecWAM | 125 | 190 |

### 5.2 Global test

The global test extends the regression case `etopo1_oper_an_fc_O48` of ecWAM to 24 h on the 36 × 29 spectral grid, from
12 UTC 31 December 2022 to 12 UTC 1 January 2023. It starts from rest with the ECMWF operational 10-m winds and sea ice, the
physics of iphys = 1 and time steps of 900 s. The errors ε_Hs and ε_Tz are the relative differences from the run with WRT,
summed over the hourly outputs from t = 1 h and over the sea points where H_s of WRT is at least 0.1 m.

| Method | ε_Hs (%) | ε_Tz (%) | Mean H_s at 24 h (m) |
|---|---|---|---|
| WRT | – | – | 1.20 |
| SWRT with K = 31 | 0.7 | 0.8 | 1.19 |
| DIA | 15 | 9 | 1.35 |

SWRT keeps the significant wave height within 1 % of WRT at every hour and nearly everywhere, while DIA gives 15 % higher
values (Figs. 1 and 2).

![Fig. 1](figs/fig_ec_global_bulk.png)

**Fig. 1.** Global statistics of the 24-h ecWAM test with S_nl from WRT (thick gray line), DIA (dashed) and SWRT with K = 31
(circles). (a) Mean H_s over the sea points; (b) ε_Hs and (c) ε_Tz at each hour.

![Fig. 2](figs/fig_ec_global_maps.png)

**Fig. 2.** ε_Hs at each sea point over the 24-h ecWAM test for (a) DIA and (b) SWRT with K = 31, with ε_Hs over all points
under each map. Seas without a value are white, mostly the ice-covered polar seas.

### 5.3 Cost

The cost is measured with DR_HOOK in the optimized build of the 24-h runs, on one of four MPI tasks, and is given relative to
DIA. WRT runs through the same vectorized routine as SWRT, and the difference in cost therefore comes from the numbers of
configurations and offsets alone.

| Method | S_nl | Time step | Run |
|---|---|---|---|
| SWRT with K = 31 | 4.0 | 1.6 | 1.4 |
| WRT | 68 | 14 | 11 |

S_nl takes about a fifth of a time step with DIA, half with SWRT and nearly all of it with WRT (Fig. 3).

![Fig. 3](figs/fig_ec_cost.png)

**Fig. 3.** Cost of S_nl in the 24-h ecWAM test. (a) Cost of SWRT with K = 31 (white) and WRT (black) relative to DIA (dashed
line) for one S_nl call, one time step and the run; (b) time per time step, split into S_nl, the other source terms and
propagation.

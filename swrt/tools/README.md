# Making the SWRT files for another spectral grid

`make_swrt_files.sh` builds the two SWRT files of section 3 for any spectral grid of ecWAM: SWRT (K = 31 configurations) and
SWRT Lite (K = 7 on at most 7 offsets). It uses training data only, the synthetic spectra of `make_synthetic.py`, and makes the
files with the steps that made the files in `swrt/files`.

## Requirements

- gfortran, and Python 3 with numpy (including f2py, with meson and ninja for its build).
- Internet access on the first run. The WRT reference is Van Vledder's Xnl code (LGPL), which is not part of this fork.
  `build_xnl.sh` downloads its four source files (`mod_constants.f90`, `mod_fileio.f90`, `serv_xnl4v5.f90`, `mod_xnl4v5.f90`)
  from WAVEWATCH III release 6.07.1 (`NOAA-EMC/WW3`, folder `model/ftn`) and builds the Python extension `xnl_wrt`.
- About 100 MB of disk in the work folder.

## Command

```bash
./make_swrt_files.sh <NANG> <NFRE> <fr1> <ifre1> [work folder]
./make_swrt_files.sh 24 36 4.177248E-02 1 ./work_d24      # the grid of case grib_restart_O96
```

NANG and NFRE are ecWAM's numbers of directions and frequencies. fr1 and ifre1 are its namelist values (FR(IFRE1) = fr1). The
frequency ratio is ecWAM's fixed 1.1. The work folder defaults to `./swrt_work_<tag>`. The tag is `d<NANG>[_n<NFRE>]_f<first
frequency>`, as `d24_f0418` or `d24_n40_f0345`. The files are written to `<work folder>/files`:

- `swrt_<tag>_K31L<L>.ww3op`: SWRT, L offsets (13 on the cost line below);
- `swrt_<tag>_K7L<L>.ww3op`: SWRT Lite;
- `choice.log`: each set tried, its training error e2 and whether it stayed finite.

`make_swrt_files.log` in the work folder gives the time of each step. Running the script again on the same work folder skips the
finished steps. It uses at most 5 processes at a time; `GRAM_JOBS` (default 4) sets the number of Gram jobs.

## Run time

On the 24 × 36 grid (d24_f0418), on a 16-core machine with another 2-process job running, the run took 16 minutes
(949 s) on 2026-11-08:

| Step | Time (s) |
|---|---|
| Xnl download and build | 5 |
| Xnl constants, WRT table | < 1 |
| Synthetic spectra (11 draws) | 15 |
| Held-out set for `check_table.py` | 6 |
| Grams of the 10 draws (4 jobs), with `check_table.py` alongside | 845 |
| SWRT: K = 31 searches on the 4 × DIA line (3 jobs) | 18 |
| SWRT: operator and divergence filter | 20 |
| SWRT Lite: search at 1.64 × DIA (4 jobs) | 35 |
| SWRT Lite: operator and divergence filter | 5 |

The Grams take most of the time. Their cost grows with NANG × NFRE.

## Selection rule

Both files are chosen from training data only. The candidates are the sets of SWRT's search on the summed Gram of the 10 synthetic
practice draws, ordered by their training error e2. SWRT takes the K = 31 sets on the line of 4 times DIA's S_nl cost
(`time_sets_ec.json`, ecWAM's cost model of the 36 × 36 grid, used for every grid). SWRT Lite takes the K = 7 sets on at most 7
offsets, searched at the cost of K = 7 on 7 offsets (1.64 times DIA). The divergence filter integrates dE/dt = S_nl(E) with
each set for 100 steps from 36 training spectra (`synfast/dyn_syn.py`). A set that does not stay finite is skipped, and the
next set in e2 order is tried. The first set that stays finite is written. When
no set of a K stays finite, the script stops with an error. The three search seeds often find the same K = 31 set, so there may
be only one candidate.

## Files

- `make_swrt_files.sh`: the entry point.
- `setup_work.py`: makes the work folder and its `grid.py` (directions, frequencies, WRT's window, JREF).
- `build_xnl.sh`: downloads the Xnl sources and builds `ext/xnl_wrt*.so` in the work folder.
- `time_sets_ec.json`: the cost model (copy of `swrt_op36/cost_ec/time_sets_ec.json`).
- `pipeline/`: the scripts copied into the work folder. They are the scripts of `swrt_grids/<tag>` and `swrt_op36/synfast`;
  each adapted script says what changed in its header.

Xnl writes WRT's quadrature file `quad99999.bqf` for the grid in the work folder at its first call (step "Xnl constants").
`check_table.py` runs as a sanity check: the full table must match fresh WRT to about 1e-4 of the largest |S_nl|.

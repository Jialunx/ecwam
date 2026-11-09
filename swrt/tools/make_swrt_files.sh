#!/bin/bash
# The two SWRT files of one spectral grid of ecWAM, from training data only, by the steps of swrt_grids/chain.sh (SWRT, K = 31) and
# swrt_op36/synfast/k_ladder.sh (SWRT Lite, K = 7 on at most 7 offsets), with the integration filter of synfast/dyn_syn.py:
#   1. work folder (setup_work.py: pipeline/, cost_ec/time_sets_ec.json, grid.py) and Xnl (build_xnl.sh, WW3 6.07.1 sources)
#   2. Xnl constants (Xnl writes quad99999.bqf for this grid), the WRT table (build_table.py), the synthetic spectra (seeds 11, 1-9, 101)
#      and the held-out set of seed 101 (testsets.py, used only by check_table.py)
#   3. check_table.py in the background (the full table against fresh WRT, a sanity check) and the Grams of the 10 practice draws
#   4. SWRT: the searches for K = 31 on the 4 x DIA line (lines_ec.py), then the distinct sets in e2 order (make_swrt_ec4x.py 31 <rank>)
#      until one stays finite in the integration on training spectra (dyn_syn.py)
#   5. SWRT Lite: search_budget.py at k_ladder.sh's budget for K = 7, L = 7, then the sets with K = 7 on at most 7 offsets in e2 order
#      (make_budget_ops.py) until one stays finite (dyn_syn.py)
#   6. files/swrt_<tag>_K31L<L>.ww3op, files/swrt_<tag>_K7L<L>.ww3op and files/choice.log
# A finished step is skipped when the script is run again on the same work folder. At most 5 processes at a time (GRAM_JOBS Grams
# and check_table.py).
# usage: ./make_swrt_files.sh <NANG> <NFRE> <fr1> <ifre1> [work folder, default ./swrt_work_<tag>]
set -o pipefail
TOOLS=$(cd "$(dirname "$0")" && pwd); [ $# -ge 4 ] || { sed -n 2,17p "$0"; exit 1; }
NA=$1; NK=$2; FR1=$3; IFRE1=$4
TAG=$(python3 -c "import sys; sys.path.insert(0, '$TOOLS'); import setup_work as s; print(s.tag_of($NA, $NK, $IFRE1, $FR1))") || exit 1
W=${5:-$(pwd)/swrt_work_$TAG}; mkdir -p "$W"; W=$(cd "$W" && pwd)
case "$W" in /tmp|/tmp/*) echo "work folder $W is under /tmp; give another folder"; exit 1;; esac
exec > >(tee -a "$W/make_swrt_files.log") 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GRAM_JOBS=${GRAM_JOBS:-4} DYN_WORKERS=1
T0=$(date +%s); TS=$T0
step () { local t=$(date +%s); echo "$(date '+%F %T') [$TAG] $* (previous step $((t - TS)) s)"; TS=$t; }
die () { step "FAILED: $*"; exit 1; }
CH=$W/files/choice.log

step "work folder $W"
python3 "$TOOLS/setup_work.py" $NA $NK $FR1 $IFRE1 "$W" > /dev/null || die "setup_work.py"
"$TOOLS/build_xnl.sh" "$W" || die "build_xnl.sh"
cd "$W"
step "Xnl constants"
[ -s quad99999.bqf ] || python3 -c "from wrt_ref import consts; consts()" >> steps.log 2>&1 || die "Xnl constants"
[ -s quad99999.bqf ] || die "Xnl wrote no quad99999.bqf"
step "table"
[ -s data/table.npz ] || python3 build_table.py >> steps.log 2>&1 || die "build_table.py"
step "synthetic spectra"
[ -s data/synthetic_s101.npz ] || python3 make_synthetic.py 11 1 2 3 4 5 6 7 8 9 101 >> steps.log 2>&1 || die "make_synthetic.py"
step "held-out set for check_table.py (WRT)"
[ -s data/test_syn_s101.npz ] || python3 testsets.py >> steps.log 2>&1 || die "testsets.py"
step "check_table.py (background) and Grams ($GRAM_JOBS jobs)"
CK=; [ -s results/check_table.json ] || { python3 check_table.py > check_table.log 2>&1 & CK=$!; }
[ -s results/gram_synthetic_s9.npz ] && [ $(ls results/gram_synthetic_s*.npz | wc -l) -eq 10 ] || ./run_grams.sh >> steps.log 2>&1 || die "Grams"
step "Grams done"
[ -n "$CK" ] && wait $CK
python3 -c "import json; d = json.load(open('results/check_table.json'))['exact']
print('check_table, full table against WRT, max |error| / max |S_nl|:', {g: f'{d[g][\"maxrel_max\"]:.1e}' for g in d})" || die "check_table.py"
step "check_table done"

[ -s $CH ] && echo "run again on this work folder, $(date '+%F %T')" >> $CH \
  || echo "SWRT files of $TAG ($NA directions x $NK frequencies, fr1 $FR1, ifre1 $IFRE1), $(date '+%F %T')" > $CH
# SWRT: K = 31 on the 4 x DIA line, the distinct sets in e2 order until one stays finite
step "searches for K = 31 on the 4x line (3 workers)"
[ -s lines_ec/line_runs_ec.jsonl ] || (cd lines_ec && python3 lines_ec.py 4 3 31 >> ../steps.log 2>&1) || die "lines_ec.py"
step "SWRT (K = 31): sets in e2 order, integration filter"
N31=; R=0
while :; do
  ls lines_ec/swrt_${TAG}_K31L*.npz > /dev/null 2>&1 && { N31=$(basename $(ls lines_ec/swrt_${TAG}_K31L*.npz) .npz); echo "K31: $N31 (finite, from an earlier run of this work folder)" | tee -a $CH; break; }
  O=$(cd lines_ec && python3 make_swrt_ec4x.py 31 $R 2>> ../steps.log) || die "make_swrt_ec4x.py 31 $R (no set left that stays finite)"
  O=$(echo "$O" | grep -m1 "^swrt_${TAG}_K31L"); N=${O%%:*}; echo "$O" >> steps.log
  python3 synfast/dyn_syn.py "$W" lines_ec/$N.npz >> steps.log 2>&1 || die "dyn_syn.py $N"
  F=$(python3 -c "import json; print(json.load(open('synfast/dyn_syn.json'))['$N']['finite'])")
  echo "K31 rank $R: $O -> $([ "$F" = True ] && echo finite || echo blows up)" | tee -a $CH
  [ "$F" = True ] && { N31=$N; break; }
  mkdir -p lines_ec/blown; for f in lines_ec/$N.npz ww3/nl6/$N.ww3op lines_ec/swrt_ec4x_K31$([ $R -gt 0 ] && echo _r$R).json; do mv $f lines_ec/blown/r${R}_$(basename $f); done
  echo "K31 rank $R blows up in the training integration; next set in e2 order" | tee -a $CH; R=$((R + 1))
done
cp ww3/nl6/$N31.ww3op files/ && echo "SWRT: $N31.ww3op" | tee -a $CH

# SWRT Lite: K = 7 on at most 7 offsets, at the budget k_ladder.sh computes for K = 7, L = 7
B=$(python3 - <<'EOF'
import json                                              # the modelled cost of K = 7 on 7 offsets, rounded up to 0.01 (k_ladder.sh)
F = json.load(open("cost_ec/time_sets_ec.json")); f = F["fit"]; L, K = 7, 7
print(f"{(f['c0'] + f['per_K'] * K + f['per_L'] * L) / F['dia_us'] + 0.005:.2f}")
EOF
)
step "SWRT Lite (K = 7, L <= 7): search at budget $B x DIA (4 workers)"
BG=$(python3 -c "print(f'{float($B):g}')"); BJ=synfast/budget_$BG.json
[ -s $BJ ] || python3 synfast/search_budget.py "$W" $B 4 > synfast/budget_$B.log 2>&1 || die "search_budget.py"
KEYS=$(python3 -c "
import json; S = json.load(open('$BJ'))['sets']; seen, out = set(), []
for k, r in sorted(S.items(), key=lambda kv: (kv[1]['e2'], kv[0])):
    u = tuple(sorted(r['units']))
    if r['K'] == 7 and r['L_limit'] <= 7 and u not in seen: seen.add(u); out.append(k)
print(' '.join(out))")
[ -n "$KEYS" ] || die "no set with K = 7 on at most 7 offsets at budget $B"
step "SWRT Lite: sets in e2 order ($KEYS), integration filter"
N7=
for k in $KEYS; do
  N=$(python3 synfast/make_budget_ops.py "$W" $B $k 2>> synfast/make_budget_ops_$B.log) || die "make_budget_ops.py $k"
  python3 -c "import json, sys; sys.exit(0 if '$N' in json.load(open('synfast/dyn_syn.json')) else 1)" 2> /dev/null \
    || python3 synfast/dyn_syn.py "$W" lines_ec/low_ops/$N.npz >> steps.log 2>&1 || die "dyn_syn.py $N"
  read F L E2 <<< $(python3 -c "import json; d = json.load(open('synfast/dyn_syn.json'))['$N']; o = json.load(open('synfast/budget_${BG}_keys.json'))['$N']
print(d['finite'], o['L'], f\"{o['e2']:.5f}\")")
  echo "K7 $k: $N, $L offsets, e2 $E2 -> $([ "$F" = True ] && echo finite || echo blows up)" | tee -a $CH
  [ "$F" = True ] && { N7=$N; break; }
  echo "K7 $k blows up in the training integration; next set in e2 order" | tee -a $CH
done
[ -n "$N7" ] || die "no set with K = 7 on at most 7 offsets stays finite"
cp ww3/nl6/$N7.ww3op files/swrt_${TAG}_K7L$L.ww3op && echo "SWRT Lite: swrt_${TAG}_K7L$L.ww3op (from $N7)" | tee -a $CH
step "done in $(( $(date +%s) - T0 )) s: $(ls files/*.ww3op | tr '\n' ' ')"

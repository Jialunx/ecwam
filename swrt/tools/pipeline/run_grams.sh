#!/bin/bash
# Grams of the 10 practice draws (build_swrt.gram, cached in results/gram_synthetic_s<seed>.npz), GRAM_JOBS (default 4) single-thread jobs
# at a time (tools copy: was 6)
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
printf "%s\n" 11 1 2 3 4 5 6 7 8 9 | xargs -P ${GRAM_JOBS:-4} -I{} python3 -c "
import os, build_swrt as b
b.gram(b.probe_set('synthetic', {}), os.path.join(b.HERE, 'results', 'gram_synthetic_s{}.npz'))
print('seed {} done', flush=True)"

#!/bin/bash
# The WRT reference of the tool: Van Vledder's Xnl (xnl4v5, LGPL) as the f2py extension xnl_wrt (modules m_constants, m_fileio,
# m_xnldata, serv_xnl4v5), built as keydoc/swrt/bench/build.sh builds xnl_bench. The four sources are not in this fork: they are
# downloaded from WAVEWATCH III release 6.07.1 (NOAA-EMC/WW3, model/ftn). f2py reads ASCII only, so the two Latin-1 degree signs in
# comments of mod_xnl4v5.f90 become "deg"; the four files are then byte-identical to the local copies of keydoc/swrt/bench.
# usage: ./build_xnl.sh <work folder>   ->  <work folder>/ext/xnl_wrt*.so (kept if it exists), sources in <work folder>/xnl_src
set -e
W=$(cd "$1" && pwd); TAG=6.07.1; URL=https://raw.githubusercontent.com/NOAA-EMC/WW3/$TAG/model/ftn
ls "$W"/ext/xnl_wrt*.so > /dev/null 2>&1 && { echo "xnl_wrt exists in $W/ext"; exit 0; }
mkdir -p "$W/xnl_src" "$W/ext"; cd "$W/xnl_src"
for f in mod_constants mod_fileio serv_xnl4v5 mod_xnl4v5; do              # module order: each file uses the ones before it
  [ -s $f.f90 ] || curl -sfL -o $f.f90 $URL/$f.f90 || { echo "download of $URL/$f.f90 failed"; exit 1; }
done
sed -i 's/\xb0/deg/g' mod_xnl4v5.f90
f2py -c mod_constants.f90 mod_fileio.f90 serv_xnl4v5.f90 mod_xnl4v5.f90 -m xnl_wrt > build_xnl.log 2>&1 || { tail -20 build_xnl.log; exit 1; }
mv xnl_wrt*.so "$W/ext/" && echo "xnl_wrt built from WW3 $TAG in $W/ext"

#!/bin/bash
# Recompute every file of groups/ from scratch (about 55 minutes in total; it can take longer on a slower computer).
# Each step prints a line with the time when it starts; the steps print nothing else while they work, so long
# silences are normal. Keep the computer awake: if it goes to sleep, the computation is paused.
# Before running it, keep the distributed data for comparison:
#     cp -r groups groups_distributed
# and afterwards compare:
#     python3 compare_outputs.py groups_distributed groups
set -e
cd "$(dirname "$0")"
if [ ! -d groups_distributed ]; then
  echo "Please make a copy of the distributed data first:  cp -r groups groups_distributed"; exit 1
fi
rm -rf groups && mkdir groups
NSTEP=13; I=0
step() { I=$((I+1)); echo; echo "=== $(date +%H:%M:%S)  [$I/$NSTEP] $*"; }
step "GAP enumeration of the groups (about 1.5 minutes)";            ./run_enum.sh
step "search and exact verification, surfaces of the paper (about 5 minutes)"
python3 batch.py groups/nec_N9.jsonl groups/nec_N12.jsonl groups/nec_N8.jsonl groups/nec_N7.jsonl \
                 groups/fuchs_N18.jsonl groups/fuchs_N12.jsonl groups/fuchs_N10.jsonl groups/fuchs_N9.jsonl \
                 groups/fuchs_N8.jsonl groups/fuchs_N7.jsonl > /dev/null
step "external test: genus 4 non-orientable, Nakamura 2009 (about 3 minutes)"
python3 batch.py groups/nec_N18.jsonl > /dev/null
step "double covers of the genus 4 surfaces, Remark 6.10 (about 10 minutes)"; python3 genus4_covers.py > /dev/null
step "tables (a few seconds)";                                    python3 summary_tables.py > /dev/null; python3 report.py > /dev/null
step "independent verification of the packings (about 1 minute; one line per family)"
python3 verify_packings.py < /dev/null
step "isometries between packings and surfaces (about 1 minute)";  python3 isometry_orbits.py > /dev/null
step "hidden automorphisms, real structures (about 7 minutes)"; python3 identify2.py > /dev/null
step "surfaces with 4 and 5 packings, X6 (about 1 minute)";        python3 fullaut.py > /dev/null; python3 cross_check_9_18.py > /dev/null
python3 weierstrass.py > /dev/null; python3 multi_extremal.py > /dev/null
step "safety margins of the floating point decisions (about 8 minutes)"; python3 margins.py > /dev/null
step "conformal automorphism groups, N=18 (about 1 minute)";  python3 conformal_aut.py fuchs_N18 > /dev/null
step "conformal automorphism groups, N=9 (about 1 minute)"; python3 conformal_aut.py fuchs_N9 > /dev/null
step "isomorphism types of the automorphism groups of all the surfaces (about 15 minutes)"; python3 aut_groups.py 2 > /dev/null
echo; echo "=== $(date +%H:%M:%S)  done. Now run:  python3 compare_outputs.py groups_distributed groups"

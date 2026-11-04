#!/bin/bash
# GAP enumeration of all the groups used in the paper (and of nec_N18, used only as an external test).
cd "$(dirname "$0")"
mkdir -p groups
run() { gap -q -b -c "N:=$1;; IDX:=$2;; ORI:=$3;; OUT:=\"groups/$4.jsonl\";;" enumerate.g >> groups/enum.log 2>&1 < /dev/null; }
run 18 18 true  fuchs_N18
run 12 24 false nec_N12
run 12 24 true  fuchs_N12
run 10 30 true  fuchs_N10
run 9  36 false nec_N9
run 9  36 true  fuchs_N9
run 8  48 false nec_N8
run 8  48 true  fuchs_N8
run 7  84 false nec_N7
run 7  84 true  fuchs_N7
run 18 36 false nec_N18        # Nakamura 2009 test (non-orientable genus 4), not part of the paper
echo ALLDONE >> groups/enum.log

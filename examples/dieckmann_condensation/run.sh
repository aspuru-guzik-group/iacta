#!/bin/sh
python rsearch.py reactants.xyz 9 32 -mtdn 120 -o claisen_cond2 -T 12 -chrg -1 -w -s 1 6 -sn 100

# Generates many many trajectories (like 1200). The product has a very low
# barrier. Low convergence rate (etemp?).

# Note also the high bond stretch.

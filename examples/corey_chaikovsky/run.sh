#!/bin/bash
python rsearch.py reactants.xyz 2 4 -mtdn 120 -o corey_chaikovsky/ -T 12 -w -s 1 6 -sn 100

# Lowest barrier belongs to the product (4 kcal /mol). Needs a higher stretch
# to get it (1x-6x).

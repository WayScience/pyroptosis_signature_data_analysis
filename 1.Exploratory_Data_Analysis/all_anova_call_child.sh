#!/bin/bash
"""
This is a script to run the anova analysis on all the cell types and features.
This is the child process that gets spun up by the parent process all_anova_call_parent.sh.
I needed to do this to get around the limit of 1000 jobs that can be submitted to the cluster.

One would call this a "cheeky" solution.
"""

feature=$1
CELL_TYPE=$2

cd scripts/

echo "Feature: $feature" "Cell type: $CELL_TYPE"

python 8.1_anova_all_groupings.py --feature $feature --cell_type $CELL_TYPE

echo "Complete"

cd ../

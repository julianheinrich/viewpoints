# This script demonstrates how to use viewpoints to create
# scenes focusing on different parts of a complex, in this case
# DNA and a p53 molecule.
#
# Assumes that viewpoints.py as been loaded such that the command
# 'set_best_view' is available.
#
# Note that you have to run all commands one-by-one instead of loading
# this script, as set_best_view returns immediately.

fetch 3igl, async = 0

as cartoon
zoom

# create selections for both p53 and DNA
select p53, chain A
select dna, chain B

# hide the selection marks
disable p53
disable dna

# set the best view for each selection and create a scene
set_best_view p53, ss, 100, 500, 500
scene p53, store

set_best_view dna, atoms, 100, 500, 500
scene dna, store

#!/usr/bin/python
__author__ = "morganlnance"

"""
HW2 Question 3 (Workshop #2 Exercise 2)

This program makes an alpha helix of a specified length
out of alanines. Phi and psi values are idealized.
This program dumps the resulting pose as 'helix.pdb'

Usage:   ./make_helix.py <number of helix residues>
Example: ./make_helix.py 20
"""


###########
# IMPORTS #
###########
import sys
from pyrosetta import init, \
    pose_from_sequence, PyMOLMover


#############
# ARGUMENTS #
#############
try:
    n_helix_residues = int(sys.argv[1])
except:
    print "\nGive me an integer for the number of helix residues.\n"
    sys.exit()


########
# MAIN #
########
# initialize pyrosetta and load the pose
init()
pose = pose_from_sequence("A"*n_helix_residues, "fa_standard")

# PyMOLMover for visualization
pmm = PyMOLMover()
pmm.keep_history(True)

# ideal phi and psi values
phi = -60
psi = -45

# set phi and psi
pmm.apply(pose)
for ii in range(1, pose.size()+1):
    pose.set_phi(ii, phi)
    pose.set_psi(ii, psi)
    pmm.apply(pose)
pose.dump_pdb("helix.pdb")

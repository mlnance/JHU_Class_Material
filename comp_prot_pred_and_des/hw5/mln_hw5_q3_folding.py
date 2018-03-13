#!/usr/bin/python
__author__="morganlnance"


"""
Morgan Nance Homework 5 Question 3 (Workshop 5 Exercise 2)

*requires 9mer and 3mer fragment files*
*file paths were hard-coded into script*

Program an ab initio folding algorithm
This program starts folding using 9mer fragments
then moves into using 3mer fragments
then alternates between Small and Shear moves
After each move, packing and minimization occurs
and the Metropolis criterion is applied
This program can send visuals to PyMOL

Usage: python <script>
Example: python mln_hw5_q3_folding.py
"""


###########
# IMPORTS #
###########
from pyrosetta import init, pose_from_sequence, get_fa_scorefxn, \
    PyMOLMover, MonteCarlo, MoveMap, standard_packer_task, \
    SequenceMover, TrialMover, RepeatMover
from pyrosetta.teaching import MinMover, PackRotamersMover
from pyrosetta import rosetta
from rosetta.core.fragment import ConstantLengthFragSet
from rosetta.protocols.simple_moves import ClassicFragmentMover, \
    SmallMover, ShearMover
from rosetta.numeric.random import random_range
from time import time


##################
# INITIALIZATION #
##################
init("-mute core.pack")

# start timing
start = time()

# RecA sequence from last 60 AA of 2reb
pose = pose_from_sequence("YKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLLSNPNSTPDFSVDDSEGVAETNEDF")

# full-atom score function
sf = get_fa_scorefxn()

# pymol mover
pmm = PyMOLMover()
pmm.keep_history(True)
pmm.apply(pose)

# MonteCarlo object for Metropolic criterion checking
kT = 1.0
mc = MonteCarlo(pose, sf, kT)

# backbone MoveMap for Small and Shear moves and fragments
mm = MoveMap()
mm.set_bb(True)

# PackRotamersMover for packing
pack_task = standard_packer_task(pose)
pack_task.restrict_to_repacking()
pack_task.or_include_current(True)
pack_mover = PackRotamersMover(sf, pack_task) 

# MinMover for minimization
min_mm = MoveMap()
min_mm.set_bb(True)
min_mm.set_chi(True)
# not setting all the other fancy options
# just focused on the MoveMap and ScoreFunction
min_mover = MinMover()
min_mover.movemap(min_mm)
min_mover.score_function(sf)


###########################
# 9MER FRAGMENT INSERTION #
###########################
print "\nFolding with 9mers...\n"
# for counting accepted moves
frag_9mer_n_accepts = 0
# set up the fragment object with the 9mers
fragset_9mer = ConstantLengthFragSet(9)
fragset_9mer.read_fragment_file("RecA_9mer_frags.frag9")
mover_9mer = ClassicFragmentMover(fragset_9mer, mm)
for ii in range(3 * pose.size()):
    # apply 9mer fragment
    mover_9mer.apply(pose)
    # pack
    pack_mover.apply(pose)
    # minimize
    min_mover.apply(pose)
    # accept or reject with Metropolic criterion
    if mc.boltzmann(pose):
        frag_9mer_n_accepts += 1
        pmm.apply(pose)


###########################
# 3MER FRAGMENT INSERTION #
###########################
print "\nFolding with 3mers...\n"
# for counting accepted moves
frag_3mer_n_accepts = 0
# set up the fragment object with the 3mers
fragset_3mer = ConstantLengthFragSet(3)
fragset_3mer.read_fragment_file("RecA_3mer_frags.frag3")
mover_3mer = ClassicFragmentMover(fragset_3mer, mm)
for ii in range(2 * pose.size()):
    # apply 3mer fragment
    mover_3mer.apply(pose)
    # pack
    pack_mover.apply(pose)
    # minimize
    min_mover.apply(pose)
    # accept or reject with Metropolic criterion
    if mc.boltzmann(pose):
        frag_3mer_n_accepts += 1
        pmm.apply(pose)


#########################
# SMALL AND SHEAR MOVES #
#########################
print "\nFolding with Small and Shear moves...\n"
n_small_and_shear_accepts = 0
# setup Small and Shear movers
small_mover = SmallMover(mm, kT, nmoves_in=5)
shear_mover = ShearMover(mm, kT, nmoves_in=5)
# either apply a small or shear move
for ii in range(2 * pose.size()):
    # True  == 1 == SmallMover
    # False == 0 == ShearMover
    if random_range(0,1):
        small_mover.apply(pose)
    else:
        shear_mover.apply(pose)
    # pack
    pack_mover.apply(pose)
    # min
    min_mover.apply(pose)
    if mc.boltzmann(pose):
        n_small_and_shear_accepts += 1
        pmm.apply(pose)

# end timing
end = time()
print "\n\nFolding simulation took %s seconds" %round((end - start),3)

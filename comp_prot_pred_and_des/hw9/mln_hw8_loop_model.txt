#!/usr/bin/python
__author__="morganlnance"


"""
Workshop 8 Question 10
Code for question 9

Couldn't really get both loops working
FoldTree wasn't doing what I expected it to
So protocol doesn't return good results
"""


from pyrosetta import init, FoldTree, MoveMap, \
    pose_from_file, get_score_function, MonteCarlo, \
    get_fa_scorefxn, PyMOLMover, PyJobDistributor
from pyrosetta import rosetta
from rosetta.protocols.loops import Loop, Loops, add_single_cutpoint_variant, \
    loop_rmsd
from rosetta.protocols.loops.loop_closure.ccd import CCDLoopClosureMover
from rosetta.core.fragment import ConstantLengthFragSet
from rosetta.protocols.simple_moves import ClassicFragmentMover
from rosetta.core.scoring import chainbreak



# initialize and load the pose
init()
pose = pose_from_file("test_in.pdb")
ref_pose = pose.clone()


# default variables
kT = 1.0
cycles = 100


# set up needed objects and movers
# ScoreFunction
sf = get_fa_scorefxn()
sf.set_weight(chainbreak, 1)
# MoveMap for the two loops of the test_in.pdb
# 15-24 (13,26,19) and 78-83 (76,85,80)
mm = MoveMap()
#[ mm.set_bb(ii, True) for ii in range(13, 26+1) ]
#[ mm.set_bb(ii, True) for ii in range(76, 85+1) ]
#[ mm.set_chi(ii, True) for ii in range(13, 26+1) ]
#[ mm.set_chi(ii, True) for ii in range(76, 85+1) ]
mm.set_bb(True)
mm.set_chi(True)
# FoldTree
ft = FoldTree()
ft.simple_tree(pose.size())
ft.new_jump(13, 26, 19)
##ft.new_jump(76, 85, 80)
pose.fold_tree(ft)
# Loops (don't include +/-2 buffer?)
loop1 = Loop(15, 25, 19)
add_single_cutpoint_variant(pose, loop1)
##loop2 = Loop(78, 83, 80)
##add_single_cutpoint_variant(pose, loop2)
##loops = Loops()
##loops.add_loop(loop1)
##loops.add_loop(loop2)
# ClassicFragmentMover
fragset_3 = ConstantLengthFragSet(3)
fragset_3.read_fragment_file("test3_fragments")
mover_3mer = ClassicFragmentMover(fragset_3, mm)
# CCD
ccd = CCDLoopClosureMover(loop1, mm)
# MonteCarlo
mc = MonteCarlo(pose, sf, kT)
# PyMOL
pmm = PyMOLMover()
pmm.keep_history(True)
pmm.apply(pose)


# extend the loops manually
loop_residues = range(15, 24)
##loop_residues.extend(range(78, 83))
for res in loop_residues:
    pose.set_phi(res, 180)
pmm.apply(pose)


# protocol
jd = PyJobDistributor("loop_output", 1, sf)
while not jd.job_complete:
    for ii in range(cycles):
        mover_3mer.apply(pose)
        ccd.apply(pose)
        pmm.apply(pose)
    loops = Loops()
    loops.add_loop(loop1)
    lrms = loop_rmsd(pose, ref_pose, loops, True)
    jd.additional_decoy_info = " LRMSD: " + str(lrms)
    jd.output_decoy(pose)

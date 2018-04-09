#!/usr/bin/python
__author__="morganlnance"

"""
Morgan Nance Homework 8 Question 2a
Workshop 4 Programming exercise 1 and 2
(but with the PyJobDistributor)
Fold 10-mer alanine or glycine chain 100 separate times
Or fold 2REB from sequence with same folding algorithm
Using the full atom folding variant with vdw and hbond scoring terms for Ala and Gly
Using the full atom score function for 2REB
Produces 100 decoys taking the lowest energy pose seen during each trajectory

Usage: python <script.py> <ALA or GLY or 2REB>
Example: python mln_hw4_q4_folding.py 2REB
"""


###########
# IMPORTS #
###########
from pyrosetta import init, PyMOLMover, \
    pose_from_sequence, get_fa_scorefxn, rosetta, Pose, \
    ScoreFunction, PyJobDistributor
from rosetta.numeric.random import gaussian, random_range, \
    uniform
from rosetta.core.scoring import hbond_sr_bb, hbond_lr_bb, \
    fa_atr, fa_rep
from rosetta.basic import periodic_range
from math import exp
import sys


#############
# ARGUMENTS #
#############
try:
    input_arg = sys.argv[1].strip().lower()
except IndexError:
    print "\nTell me what sequence you want to fold."
    print "Try typing 'gly', 'ala', or '2reb'\n"
    sys.exit()
# ensure a valid argument was passed
if input_arg not in ["gly", "ala", "2reb"]:
    print "\nI can only do GLY, ALA, or 2REB as an argument.\n"
    sys.exit()


##################
# INITIALIZATION #
##################
init("-multiple_processes_writing_to_one_directory")
if input_arg == "ala":
    pose_fresh = pose_from_sequence("A" * 10)
elif input_arg == "gly":
    pose_fresh = pose_from_sequence("G" * 10)
else:
    # sequence of 2REB from PDB
    pose_fresh = pose_from_sequence("AIDENKQKALAAALGQIEKQFGKGSIMRLGEDRSMDVETISTGSLSLDIALGAGGLPMGRIVEIYGPESSGKTTLTLQVIAAAQREGKTCAFIDAEHALDPIYARKLGVDIDNLLCSQPDTGEQALEICDALARSGAVDVIVVDSVAALTPKAEIEGEIGDSHMGLAARMMSQAMRKLAGNLKQSNTLLIFINQIRMKIGVMFGNPETTTGGNALKFYASVRLDIRRIGAVKEGENVVGSETRVKVVKNKIAAPFKQAEFQILYGEGINFYGELVDLGVKEKLIEKAGAWYSYKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLLSNPNSTPDFSVDDSEGVAETNEDF")

fa_sf = get_fa_scorefxn()
sf = ScoreFunction()
pmm = PyMOLMover()
pmm.keep_history(True)

# declarations from workshop and homework
angle_deviation = 25.0
kT = 1.0
nstruct = 10
n_cycles = 1000

# if doing ala or gly, create the appropriate sf
if input_arg in ["gly", "ala"]:
    # set the scorefunction according to the workshop
    # only van der Waals (fa_atr and fa_rep) and
    # hbonds should be used. Ala and Gly make no
    # side chain hbonds, so use the bb terms
    # get the weights from the full atom score function
    sf.set_weight(fa_atr, fa_sf.get_weight(fa_atr))
    sf.set_weight(fa_rep, fa_sf.get_weight(fa_rep))
    sf.set_weight(hbond_sr_bb, fa_sf.get_weight(hbond_sr_bb))
    sf.set_weight(hbond_lr_bb, fa_sf.get_weight(hbond_lr_bb))
# otherwise, use the standard full atom sf
else:
    sf.assign( fa_sf.clone() )


#############
# FUNCTIONS #
#############
def make_move(pose):
    """
    Make a random phi or psi move on a random residue in the pose
    :param pose: Pose
    :return: Pose
    """
    # pick a random residue from 1 to N
    rand_res = random_range(1, pose.size())
    # determine the amount of deviation we will make
    deviation = gaussian() * angle_deviation
    # make a phi or psi move
    # phi = 1 = True and psi = 0 = False
    if random_range(0, 1):
        old_phi = pose.phi(rand_res)
        # torsions are -180 to 180 so use 360 periodic range
        new_phi = periodic_range( old_phi - deviation, 360 )
        pose.set_phi(rand_res, new_phi)
    else:
        old_psi = pose.psi(rand_res)
        # torsions are -180 to 180 so use 360 periodic range
        new_psi = periodic_range( old_psi - deviation, 360 )
        pose.set_psi(rand_res, new_psi)

    return pose


def accept_move(pose_current, pose_new):
    """
    Give a pose_current, compare its score against pose_new and
    return True if the move is to be accepted
    :param pose_current: Pose
    :param pose_new: Pose
    :return: bool
    """
    # score the poses
    current_score = sf(pose_current)
    new_score = sf(pose_new)

    # compare the scores
    # if new_score < current_score, return pose_new
    if new_score < current_score:
        return True
    # otherwise, apply the Metropolis criterion
    # P = exp( -dE/kT )
    else:
        # calculate dE (new - old) which is a positive number
        dE = new_score - current_score
        # calculate the probability of accepting this pose
        # this should return a small number less than 1
        # the more positive dE is, the smaller P is
        P = exp( -dE / kT )
        # generate a random uniform number between [0, 1)
        R = uniform()
        # if R < P, accept the move and return pose_new
        # the larger dE is, the smaller P is, the smaller
        # chance that R is smaller than P, thus smaller
        # chance of accepting that move
        if R < P:
            return True
        # otherwise exit this statement and return pose_current
        # the move was not accepted

    return False


################
# FOLDING MAIN #
################
# visualize before folding
#pmm.apply(ala_pose)


# store the lowest energy pose seen for each independent trajectory
lowE_poses_seen = []

# wrap the main protocol within the PyJobDistributor
jd = PyJobDistributor("hw8_output", nstruct, sf)
jd.native_pose = pose_fresh.clone()
# while the job is not complete, keep producing decoys
# this outer while loop is for decoy trajectories
while not jd.job_complete:
    print "\n\nFolding, please wait..."
    # get a fresh copy of the pose
    pose = pose_fresh.clone()
    pose.pdb_info().name(jd.current_name)
    # for keeping track of the lowest-E pose seen each trajectory
    lowE_pose = pose_fresh.clone()

    # fold
    # this inner loop is how many moves get made per decoy
    for jj in range(n_cycles):
        # keep a copy of the pose before the move was made
        # needed for Metropolis comparison
        pose_before_move = pose.clone()

        # perturb the pose with a phi or psi move
        make_move(pose)

        # compare the pose before the move and of after
        # Metropolis function will return the better of the two poses
        # so assign the current pose to whichever one is best
        if accept_move( pose_before_move, pose ):
            # move was accepted, so visualize
            #pmm.apply(pose)
            pass
        else:
            # move was rejected, revert pose to before move
            pose.assign( pose_before_move )
            # update the lowE_pose if needed
            # this is the lowest E pose seen during folding
        if sf(pose) < sf(lowE_pose):
            lowE_pose.assign(pose)

    # end trajectory for this decoy, store lowE_pose
    lowE_poses_seen.append(lowE_pose)
    # and output the pose
    jd.output_decoy(lowE_pose)

#!/usr/bin/python
__author__="morganlnance"

"""
Morgan Nance Homework 4 Question 4
Workshop 4 Programming exercise 1
Fold 10-mer alanine chain 100 separate times
Using the full atom folding variant with vdw and hbond scoring terms

Usage: python mln_hw4_q4_folding.py
"""


###########
# IMPORTS #
###########
from pyrosetta import init, PyMOLMover, \
    pose_from_sequence, get_fa_scorefxn, rosetta, Pose, \
    ScoreFunction
from rosetta.numeric.random import gaussian, random_range, \
    uniform
from rosetta.core.scoring import hbond_sr_bb, hbond_lr_bb, \
    fa_atr, fa_rep
from rosetta.basic import periodic_range
from math import exp


##################
# INITIALIZATION #
##################
init()
ala_pose_fresh = pose_from_sequence("A" * 10)
gly_pose_fresh = pose_from_sequence("G" * 10)
fa_sf = get_fa_scorefxn()
sf = ScoreFunction()
pmm = PyMOLMover()
pmm.keep_history(True)

# declarations from workshop and homework
angle_deviation = 25.0
kT = 1.0
n_trajectories = 100
n_cycles = 100

# set the scorefunction according to the workshop
# only van der Waals (fa_atr and fa_rep) and
# hbonds should be used. Ala and Gly make no
# side chain hbonds, so use the bb terms
# get the weights from the full atom score function
sf.set_weight(fa_atr, fa_sf.get_weight(fa_atr))
sf.set_weight(fa_rep, fa_sf.get_weight(fa_rep))
sf.set_weight(hbond_sr_bb, fa_sf.get_weight(hbond_sr_bb))
sf.set_weight(hbond_lr_bb, fa_sf.get_weight(hbond_lr_bb))


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
lowE_ala_poses_seen = []
lowE_gly_poses_seen = []

for ii in range(n_trajectories):
    # get a fresh copy of the pose
    ala_pose = ala_pose_fresh.clone()
    gly_pose = gly_pose_fresh.clone()
    # for keeping track of the lowest-E pose seen each trajectory
    lowE_ala_pose = ala_pose_fresh.clone()
    lowE_gly_pose = gly_pose_fresh.clone()

    # fold
    for jj in range(n_cycles):
        # keep a copy of the pose before the move was made
        # needed for Metropolis comparison
        ala_pose_before_move = ala_pose.clone()
        gly_pose_before_move = gly_pose.clone()

        # perturb the pose with a phi or psi move
        make_move(ala_pose)
        make_move(gly_pose)

        # Alanine 10-mer
        # compare the pose before the move and of after
        # Metropolis function will return the better of the two poses
        # so assign the current pose to whichever one is best
        if accept_move( ala_pose_before_move, ala_pose ):
            # move was accepted, so visualize
            #pmm.apply(ala_pose)
            pass
        else:
            # move was rejected, revert pose to before move
            ala_pose.assign( ala_pose_before_move )
        # update the lowE_ala_pose if needed
        # this is the lowest E pose seen during folding
        if sf(ala_pose) < sf(lowE_ala_pose):
            lowE_ala_pose.assign(ala_pose)

        # Glycine 10-mer
        if accept_move( gly_pose_before_move, gly_pose ):
            #pmm.apply(gly_pose)
            pass
        else:
            gly_pose.assign( gly_pose_before_move )
        if sf(gly_pose) < sf(lowE_gly_pose):
            lowE_gly_pose.assign(gly_pose)

    # end trajectory, store lowE_ala_pose
    lowE_ala_poses_seen.append(lowE_ala_pose)
    # end trajectory, store lowE_gly_pose
    lowE_gly_poses_seen.append(lowE_gly_pose)

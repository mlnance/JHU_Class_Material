#!/usr/bin/python
__author__ = "morganlnance"


'''
'''


###########
# IMPORTS #
###########
import sys
from math import sqrt
from dihedral_angle import calc_dihedral_angle


###############
# DEFINITIONS #
###############
sugar_residue_names = ["BGC", "GLC", "BMA",
                       "MAN", "GAL", "NAG"]
# pyranose rings have 5 carbons in the ring
sugar_ring_atoms = ["C1", "C2", "C3", "C4", "C5"]
# for pruning which TRP to sugar CH atom contacts to keep
# based on distance
trp_to_sugar_ch_cutoff = 6


##############
# INPUT ARGS #
##############
try:
    in_file = sys.argv[1]
except IndexError:
    print "\n\nPlease give me a .pdb file\n"
    sys.exit()


#############
# READ FILE #
#############
with open(in_file, "r") as fh:
    pdb_lines = fh.readlines()
# pull out desired lines
atom_lines = [l for l in pdb_lines
              if l.startswith("ATOM")]
hetatm_lines = [l for l in pdb_lines
                if l.startswith("HETATM")
                and "HOH" not in l]
trp_lines = [l for l in atom_lines
             if "TRP" in l]
# delete the pdb_lines now
del(pdb_lines)


##################################
# BUILD A TRP RESIDUE DICTIONARY #
##################################
# key: resnum_chain, value: {x, y, z}
trp_cd2_ce2_midpoints = {}
# build a trp dictionary for each residue's lines
# key: resnum_chain, value: [lines]
trp_dict = {}
for trp_line in trp_lines:
    # build the id for the trp residue
    resnum = trp_line[23:26].strip()
    chain = trp_line[21:22]
    trp_id = resnum + "_" + chain
    # if this id isn't in the dictionary
    # add the id as the key and the line to a list
    if trp_id not in trp_dict.keys():
        trp_dict[trp_id] = [trp_line]
    # if we already have this trp_id, add the trp_line
    else:
        trp_dict[trp_id].append(trp_line)        


#############################
# GET TRP CD2 CE2 MIDPOINTS #
#############################
# key: resnum_chain, value: [x, y, z]
trp_cd2_ce2_midpoints = {}
# trp_id = resnum_chain
for trp_id, trp_res in trp_dict.iteritems():
    cd2_line = ""
    ce2_line = ""
    # get the CD2 and CE2 atom info
    for atom_line in trp_res:
        if atom_line[12:16].strip() == "CD2":
            # xyz info for CD2 atom
            cd2_line = atom_line
            cd2_x = float(atom_line[30:38].strip())
            cd2_y = float(atom_line[38:46].strip())
            cd2_z = float(atom_line[46:54].strip())
        if atom_line[12:16].strip() == "CE2":
            # xyz info for CE2 atom
            ce2_line = atom_line
            ce2_x = float(atom_line[30:38].strip())
            ce2_y = float(atom_line[38:46].strip())
            ce2_z = float(atom_line[46:54].strip())
    # calculate and round the midpoint between CD2 and CE2
    midpoint = [ round((cd2_x + ce2_x)/2, 3), 
                 round((cd2_y + ce2_y)/2, 3),
                 round((cd2_z + ce2_z)/2, 3) ]
    # add to the dictionary
    trp_cd2_ce2_midpoints[trp_id] = midpoint


####################################
# BUILD A SUGAR RESIDUE DICTIONARY #
####################################
# key: resnum_chain, value: {x, y, z}
sugar_cd2_ce2_midpoints = {}
# build a sugar dictionary for each residue's lines
# key: resnum_chain, value: [lines]
sugar_dict = {}
for sugar_line in hetatm_lines:
    # check that this residue is in the defined set of sugar residues
    resname = sugar_line[17:20].strip()
    if resname not in sugar_residue_names:
        pass
    # build the id for the sugar residue
    resnum = sugar_line[23:26].strip()
    chain = sugar_line[21:22]
    sugar_id = resnum + "_" + chain
    # if this id isn't in the dictionary
    # add the id as the key and the line to a list
    if sugar_id not in sugar_dict.keys():
        sugar_dict[sugar_id] = [sugar_line]
    # if we already have this sugar_id, add the sugar_line
    else:
        sugar_dict[sugar_id].append(sugar_line)


###############################################
# FIND CLOSEST SUGAR CH ATOMS TO TRP MIDPOINT #
###############################################
# trp to sugar dictionary
# key: TRPresnum_chain+SUGARresnum_chain_atomname, value: distance
# example: 127_A+310_A_C3 : 5.2
all_trp_midpt_to_sugar_ch_dict = {}
# for each tryptophan, calc distances to each sugar CH
# then determine which one is the minimum
for trp_id, trp_midpt in trp_cd2_ce2_midpoints.iteritems():
    midpt_to_ch_atom_dists = {}
    # pull out the TRP midpoint x, y, z
    trp_x = trp_cd2_ce2_midpoints[trp_id][0]
    trp_y = trp_cd2_ce2_midpoints[trp_id][1]
    trp_z = trp_cd2_ce2_midpoints[trp_id][2]
    # iterate through each sugar residue
    # and each CH atom in that sugar
    for sugar_id, sugar_res in sugar_dict.iteritems():
        for atom_line in sugar_res:
            # check only the relevant CH atoms
            atom_name = atom_line[12:16].strip()
            if atom_name in sugar_ring_atoms:
                # add onto the sugar_id to make sugar_atom_id
                # resnum_chain_atomname
                sugar_atom_id = sugar_id + "_" + atom_name
                # calculate the distance to the TRP midpoint
                ch_x = float(atom_line[30:38].strip())
                ch_y = float(atom_line[38:46].strip())
                ch_z = float(atom_line[46:54].strip())
                # distance calc
                # sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
                dist = sqrt( (ch_x - trp_x)**2 +
                             (ch_y - trp_y)**2 +
                             (ch_z - trp_z)**2 )
                # add the midpoint to CH distance to the list
                # for this current TRP's midpoint
                midpt_to_ch_atom_dists[sugar_atom_id] = dist

    # done checking all sugar CH atoms for this specific midpoint
    # determine which distance is the minimum
    # and update the data dictionary accordingly
    min_dist = None
    min_sugar_atom_id = None
    for sugar_atom_id, dist in midpt_to_ch_atom_dists.iteritems():
        # if this is the first distance checked
        # or if this distance is smaller than it if it's not None
        # update the holders
        if (min_dist is None) or (dist < min_dist):
            min_dist = dist
            min_sugar_atom_id = sugar_atom_id

    # add this data to the dictionary
    # key: TRPresnum_chain+SUGARresnum_chain_atomname, value: distance
    # example: 127_A+310_A_C3 : 5.2
    trp_to_sugar_atom_id = trp_id + "+" + min_sugar_atom_id
    all_trp_midpt_to_sugar_ch_dict[trp_to_sugar_atom_id] = min_dist

# we collected all the minimum TRP midpoints to sugar CH atom
# ie) each TRP in the prot must have some minimum distance
# to some sugar CH atom in the protein. that distance could be
# large though. So go through the dictionary and prune out
# TRP to sugar CH distances below some cutoff distance
# (iterate through the same dictionary key and value pairs
# and keep only values that are below the cutoff)
# *** these are the TRP to sugar CH atom contacts we want to futher analyze
trp_midpt_to_sugar_ch_dict = { trp_to_sugar_atom_id:
                                   round(min_dist, 3)
                               for trp_to_sugar_atom_id, min_dist 
                               in all_trp_midpt_to_sugar_ch_dict.iteritems()
                               if min_dist < trp_to_sugar_ch_cutoff }


#############################################
# CALCULATE DIHEDRAL BETWEEN NE1 CE2 CD2 CH #
#############################################

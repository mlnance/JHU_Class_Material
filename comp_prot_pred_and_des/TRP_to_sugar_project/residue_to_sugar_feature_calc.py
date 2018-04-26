#!/usr/bin/python


'''
TODO: make this work on a directory of pdb structs and write to one df
have to change how data gets added to df? dictionaries are instantiated when needed
maybe make all dicts in the beginning before for loop and then add to them continuously
store pdb name in the df too
'''


###########
# IMPORTS #
###########
import sys
from math import sqrt
import pandas as pd
from dihedral_angle import calc_dihedral_angle


###############
# DEFINITIONS #
###############
sugar_residue_names = ["BGC", "GLC", "BMA",
                       "MAN", "GAL", "NAG"]
# pyranose rings have 5 carbons in the ring
sugar_ring_atoms = ["C1", "C2", "C3", "C4", "C5"]
# for pruning which TRP to sugar CH atom contacts to keep
# based on distance in Angstroms
trp_to_sugar_ch_cutoff = 6


##############
# INPUT ARGS #
##############
try:
    in_file = sys.argv[1]
except IndexError:
    print "\n\nPlease give me a .pdb file\n"
    sys.exit()


'''
if the input is a pdb file, put it into a list
if the input is a directory, put all .pdb files into a list
have to then put everything below into a for loop, requiring one tab over
'''


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
# need to also keep track of the xyz coordinates of the CH atom
# that is the closest to the TRP midpoint atom
# key: SUGARresnum_chain_atomname, value: [x, y, z] of CH atom
all_sugar_ch_atom_xyz_dict = {}
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
                # add the sugar CH atom xyz coords to the dict
                all_sugar_ch_atom_xyz_dict[sugar_atom_id] = [ch_x,
                                                             ch_y,
                                                             ch_z]
                

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
# ie) each TRP in the protein must have some minimum distance
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
# with the pruned TRP midpoint data, we can figure out which CH atom was closest
# and prune a dictionary of sugar CH atoms to their xyz coords for later use
# trp_to_sugar_atom_id looks like: TRPresnum_chain+SUGARresnum_chain_atomname
# so split on the "+" and pull out the sugar id to get its xyz coords
sugar_ch_atom_xyz_dict = { trp_to_sugar_atom_id.split("+")[-1]:
                               all_sugar_ch_atom_xyz_dict[trp_to_sugar_atom_id.split("+")[-1]]
                           for trp_to_sugar_atom_id
                           in trp_midpt_to_sugar_ch_dict.keys() }


#######################################################
# CALCULATE DIHEDRAL BETWEEN TRP NE1 CE2 CD2 SUGAR CH #
#######################################################
# for each TRP residue that has a relevant dihedral to calculate
# pull out the TRP residue from the trp_dict, get the relevant atoms
# and then pull out the corresponding sugar CH atom from the
# sugar_ch_atom_xyz_dict dictionary
# key: TRPresnum_chain+SUGARresnum_chain_atomname, value: NE1-CE2-CD2-CH dihedral
trp_to_sugar_ch_dihedral_dict = {}
for trp_to_sugar_atom_id in trp_midpt_to_sugar_ch_dict.keys():
    # pull out the TRPresnum_chain id
    trp_resnum_chain = trp_to_sugar_atom_id.split("+")[0]
    # pull out the sugar CH atom id too
    sugar_atom_id = trp_to_sugar_atom_id.split("+")[-1]

    # grab the TRP residue lines given the trp_resnum_chain
    trp_res_lines = trp_dict[trp_resnum_chain]
    # pull out the NE1 CE2 and CD2 atom coordinates
    # start them out as None in case for some reason the TRP
    # residue doesn't have these atoms and we can't calculate the dihedral
    trp_NE1_xyz = None
    trp_CE2_xyz = None
    trp_CD2_xyz = None
    for trp_line in trp_res_lines:
        if trp_line[12:16].strip() == "NE1":
            trp_NE1_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]
        if trp_line[12:16].strip() == "CE2":
            trp_CE2_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]
        if trp_line[12:16].strip() == "CD2":
            trp_CD2_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]

    # grab the xyz coordinates of the relevant sugar CH atom
    # that is closest to this TRP's midpt atom
    sugar_ch_atom_xyz = sugar_ch_atom_xyz_dict[sugar_atom_id]

    # calculate the dihedral between NE1 - CE2 - CD2 - sugar CH atom
    trp_to_sugar_dihedral = round(calc_dihedral_angle(trp_NE1_xyz,
                                                      trp_CE2_xyz,
                                                      trp_CD2_xyz,
                                                      sugar_ch_atom_xyz), 3)
    # add to the dictionary given the TRP to sugar CH id key
    trp_to_sugar_ch_dihedral_dict[trp_to_sugar_atom_id] = trp_to_sugar_dihedral


###############################################
# CALCULATE DIHEDRAL BETWEEN TRP CA CB CG CD1 #
###############################################
# for all TRP residues that are CH-pi interacting based on the above criteria
# calculate its chi2 value from atoms CA CB CG and CD1
trp_chi2_dihedral_dict = {}
for trp_to_sugar_atom_id in trp_midpt_to_sugar_ch_dict.keys():
    # pull out the TRPresnum_chain id
    trp_resnum_chain = trp_to_sugar_atom_id.split("+")[0]

    # grab the TRP residue lines given the trp_resnum_chain
    trp_res_lines = trp_dict[trp_resnum_chain]
    # pull out the CA CB CG and CD1 atom coordinates
    # start them out as None in case for some reason the TRP
    # residue doesn't have these atoms and we can't calculate the dihedral
    trp_CA_xyz = None
    trp_CB_xyz = None
    trp_CG_xyz = None
    trp_CD1_xyz = None
    for trp_line in trp_res_lines:
        if trp_line[12:16].strip() == "CA":
            trp_CA_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]
        if trp_line[12:16].strip() == "CB":
            trp_CB_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]
        if trp_line[12:16].strip() == "CG":
            trp_CG_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]
        if trp_line[12:16].strip() == "CD1":
            trp_CD1_xyz = [float(trp_line[30:38].strip()),
                           float(trp_line[38:46].strip()),
                           float(trp_line[46:54].strip())]

    # calculate the dihedral between CA - CB - CG - CD1
    trp_chi2_dihedral = round(calc_dihedral_angle(trp_CA_xyz,
                                                  trp_CB_xyz,
                                                  trp_CG_xyz,
                                                  trp_CD1_xyz), 3)
    # add to the dictionary given the TRPresnum_chain
    trp_chi2_dihedral_dict[trp_resnum_chain] = trp_chi2_dihedral


###################
# DATA COLLECTION #
###################
# build the dataframe column by column
df = pd.DataFrame()
# TRP residue number and chain
# pull from dictionary that has TRP and sugar id together
# to ensure the order of collecting TRP and sugar ids corresponds
# to the order that the TRP and sugar residues are actually interacting
# key: TRPresnum_chain+SUGARresnum_chain_atomname
trp_to_sugar_ids = trp_midpt_to_sugar_ch_dict.keys()
trp_ids = [trp_to_sugar_id.split("+")[0] for trp_to_sugar_id in trp_to_sugar_ids]
trp_resnums = [trp_id.split("_")[0] for trp_id in trp_ids]
trp_chains = [trp_id.split("_")[1] for trp_id in trp_ids]
# sugar CH residue number, chain, and CH atom name
# sugar_id: SUGARresnum_chain_atomname
sugar_ids = [trp_to_sugar_id.split("+")[1] for trp_to_sugar_id in trp_to_sugar_ids]
sugar_resnums = [sugar_id.split("_")[0] for sugar_id in sugar_ids]
sugar_chains = [sugar_id.split("_")[1] for sugar_id in sugar_ids]
sugar_ch_atoms = [sugar_id.split("_")[2] for sugar_id in sugar_ids]
# TRP midpoint [x, y, z] and distances to sugar CH atom
trp_midpt_xyzs = [trp_cd2_ce2_midpoints[trp_id] for trp_id in trp_ids]
# a little redundant, but explicit to ensure the correct order of things
trp_midpt_to_ch_dists = [trp_midpt_to_sugar_ch_dict[trp_to_sugar_id]
                         for trp_to_sugar_id in trp_to_sugar_ids]
# TRP NE1 CE2 CD2 to sugar CH atom dihedrals
# also a little redundant, but again, ensuring correct order
trp_to_sugar_ch_dihedrals = [trp_to_sugar_ch_dihedral_dict[trp_to_sugar_id]
                             for trp_to_sugar_id in trp_to_sugar_ids]
# TRP chi2 dihedrals
trp_chi2_dihedrals = [trp_chi2_dihedral_dict[trp_id] for trp_id in trp_ids]


# add data to the dataframe
df["TRP_resnum"] = trp_resnums
df["TRP_chain"] = trp_chains
df["sugar_resnum"] = sugar_resnums
df["sugar_chain"] = sugar_chains
df["sugar_ch_atom"] = sugar_ch_atoms
df["TRP_midpt_xyz"] = trp_midpt_xyzs
df["TRP_midpt_to_CH_dist"] = trp_midpt_to_ch_dists
df["TRP_to_CH_dihedral"] = trp_to_sugar_ch_dihedrals
df["TRP_chi2_dihedral"] = trp_chi2_dihedrals

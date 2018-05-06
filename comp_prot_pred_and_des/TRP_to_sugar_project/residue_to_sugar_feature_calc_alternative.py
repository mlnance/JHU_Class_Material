#!/usr/bin/python


'''
TODO: make this work on a directory of pdb structs and write to one df
have to change how data gets added to df? dictionaries are instantiated when needed
maybe make all dicts in the beginning before for loop and then add to them continuously
store pdb name in the df too
clarify intent of getting atom names and stuff
clarify why I am using an atom id
edit so that no multi-occupancy counting occurs (ensure blank or A)
three TRP atoms to 1 sugar atom doesn't really capture planarity though
look at 1tfv as an example
how to check planarity between trp and sugar?
'''


###########
# IMPORTS #
###########
import sys
import os
import gzip
from math import sqrt
from math import pi
from time import time
import pandas as pd
from dihedral_angle import calc_dihedral_angle
from sympy import Plane, Point3D


###############
# DEFINITIONS #
###############
sugar_residue_names = ["BGC", "GLC", "BMA",
                       "MAN", "GAL", "NAG",
                       "FUC",
                       "Glc", "Man", "Fuc"]
# pyranose rings have 5 carbons in the ring
# and one oxygen atom
# I'm sure there are exceptions, but first round should be generic
# and specific to get an understanding of the data
sugar_ring_atoms = ["C1", "C2", "C3", "C4", "C5", "O5"]
# for pruning which TRP to sugar midpoint contacts to keep
# based on distance in Angstroms
trp_to_sugar_midpt_cutoff = 4.5
# dataframe for data
df_columns = ["PDB", "TRP_resnum", "TRP_chain",
              "sugar_resname", "sugar_resnum", "sugar_chain",
              "TRP_midpt_xyz", "sugar_midpt_xyz",
              "TRP_midpt_to_sugar_midpt_dist",
              "TRP_to_sugar_plane_angle_raw",
              "TRP_to_sugar_plane_angle",
              "TRP_chi2_dihedral"]
df = pd.DataFrame(columns = df_columns)
# start time
start = time()


##############
# INPUT ARGS #
##############
try:
    user_input = sys.argv[1]
except IndexError:
    print "\n\nPlease give me a .pdb(.gz) file or a directory of files\n"
    sys.exit()
# is the input a single .pdb (or .pdb.gz) file?
if os.path.isfile(user_input) and ".pdb" in user_input:
    # put the file into a single list
    pdb_files = [os.path.abspath(user_input)]
elif os.path.isdir(user_input):
    # put the absolute path of all .pdb(.gz) into a list
    pdb_files = [os.path.abspath(user_input) + "/" + f
                 for f in os.listdir(user_input)
                 if ".pdb" in f]
else:
    print "\n\nPlease give me a .pdb(.gz) file or a directory of files\n"
    sys.exit()


######################
# OPEN AND READ FILE #
######################
# for every pdb file
for pdb_file in pdb_files:
    # open the file and get the lines
    # if it's a .pdb file, open normally
    if pdb_file.endswith(".pdb"):
        try:
            with open(pdb_file, "r") as fh:
                pdb_lines = fh.readlines()
        except:
            # if this file doesn't open and read, just skip it
            continue
    # if it's a .pdb.gz file, open with gzip.open
    elif pdb_file.endswith(".pdb.gz"):
        try:
            with gzip.open(pdb_file, "r") as fh:
                pdb_lines = fh.readlines()
        except:
            # if this file doesn't open and read, just skip it
            continue
    # else, not sure how this file got here, continue
    else:
        continue
    # get the name of the PDB file
    # example: /path/to/file/1abc_long_rosetta_suffix.pdb.gz ---> 1abc
    pdb_id = pdb_file.split("/")[-1].split(".")[0].split("_")[0]
    # example: /path/to/file/1abc_suffix.pdb ---> 1abc_suffix
    pdb_name = pdb_file.split("/")[-1].split(".")[0]

    # pull out desired ATOM and HETATM lines
    # ***** IMPORTANT NOTE *****
    # if it's a regular PDB file, sugar residues will be a HETATM
    # if it's a Rosetta PDB file, sugar residues will be an ATOM
    atom_lines = [l for l in pdb_lines
                  if l.startswith("ATOM")
                  or l.startswith("HETATM") and l[17:20].strip() != "HOH"]
    # pull out just the TRP residue lines for quicker parsing later
    # and only pull out TRP residues of one conformation
    # so either they are full occupancy or take the A occupancy
    trp_lines = [l for l in atom_lines
                 if l[17:20].strip() == "TRP"
                 and (l[16:17] == " " or l[16:17] == "A")]
    # delete the pdb_lines now as we don't need all of them
    # and might as well free up some space
    del(pdb_lines)


##################################
# BUILD A TRP RESIDUE DICTIONARY #
##################################
    # build a trp dictionary for each residue's lines
    # key: PDBname_TRPresnum_chain, value: [lines]
    trp_dict = {}
    for trp_line in trp_lines:
        # build the id for the TRP residue
        # each unique TRP should have a unique TRP id
        resnum = trp_line[22:26].strip()
        chain = trp_line[21:22]
        trp_id = pdb_id + "_" + resnum + "_" + chain
        # if this id isn't in the dictionary
        # add the id as the key and the line to a list
        if trp_id not in trp_dict.keys():
            trp_dict[trp_id] = [trp_line]
        # if we already have this trp_id, add the trp_line
        else:
            trp_dict[trp_id].append(trp_line)
    # if no TRP residues, skip this PDB
    if len(trp_dict) == 0:
        continue

#############################
# GET TRP CD2 CE2 MIDPOINTS #
#############################
    # key: PDBname_TRPresnum_chain, value: [x, y, z]
    trp_cd2_ce2_midpoints = {}
    # trp_id = PDBname_TRPresnum_chain
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
    # build a sugar dictionary for each residue's lines
    # key: PDBname_SUGARresname_resnum_chain, value: [lines]
    sugar_dict = {}
    for sugar_line in atom_lines:
        # check that this residue is in the defined set of sugar residues
        resname = sugar_line[17:20].strip()
        if resname not in sugar_residue_names:
            continue
        # build the id for the sugar residue
        resnum = sugar_line[22:26].strip()
        chain = sugar_line[21:22]
        sugar_id = pdb_id + "_" + resname + "_" + resnum + "_" + chain
        # if this id isn't in the dictionary
        # add the id as the key and the line to a list
        if sugar_id not in sugar_dict.keys():
            sugar_dict[sugar_id] = [sugar_line]
        # if we already have this sugar_id, add the sugar_line
        else:
            sugar_dict[sugar_id].append(sugar_line)
    # if no sugar residues, skip this PDB
    if len(sugar_dict) == 0:
        continue


#################################
# GET SUGAR O5 AND C3 MIDPOINTS #
#################################
    # key: PDBname_SUGARresname_resnum_chain, value: [x, y, z]
    sugar_o5_c3_midpoints = {}
    # sugar_id = PDBname_SUGARresname_resnum_chain
    for sugar_id, sugar_res in sugar_dict.iteritems():
        o5_line = ""
        c3_line = ""
        # get the O5 and C3 atom info
        for atom_line in sugar_res:
            if atom_line[12:16].strip() == "O5":
                # xyz info for O5 atom
                o5_line = atom_line
                o5_x = float(atom_line[30:38].strip())
                o5_y = float(atom_line[38:46].strip())
                o5_z = float(atom_line[46:54].strip())
            if atom_line[12:16].strip() == "C3":
                # xyz info for C3 atom
                c3_line = atom_line
                c3_x = float(atom_line[30:38].strip())
                c3_y = float(atom_line[38:46].strip())
                c3_z = float(atom_line[46:54].strip())
        # calculate and round the midpoint between O5 and C3
        midpoint = [ round((o5_x + c3_x)/2, 3), 
                     round((o5_y + c3_y)/2, 3),
                     round((o5_z + c3_z)/2, 3) ]
        # add to the dictionary
        sugar_o5_c3_midpoints[sugar_id] = midpoint


#################################################
# FIND DISTANCE BETWEEN TRP AND SUGAR MIDPOINTS #
#################################################
    # TRP midpoint to sugar midpoint atom distance dictionary
    # key: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain
    # value: distance
    # example: 3ACH_127_A+3ACH_Glc_310_A : 5.2
    all_trp_midpt_to_sugar_midpt_dict = {}
    # need to also keep track of the xyz coordinates of the sugar midpoints
    # that is the closest to the TRP midpoint atom
    # key: PDBname_SUGARresname_resnum_chain, value: [x, y, z] of the midpoint
    # collecting information for ALL TRP to sugar midpoint calculations
    # then will prune after collecting this data to only keep values
    # that are within some relevant cutoff range
    all_sugar_midpt_xyz_dict = {}
    # for each tryptophan midpoint, calc distances to each sugar midpoint
    # then determine which one is the minimum
    for trp_id, trp_midpt in trp_cd2_ce2_midpoints.iteritems():
        # each unique TRP will have a dictionary of its midpoint
        # to all sugar midpoints. that will then be pruned
        # for contacts that are within some relevant cutoff range
        trp_midpt_to_sugar_midpt_dists = {}
        # pull out the TRP midpoint x, y, z
        trp_x = trp_cd2_ce2_midpoints[trp_id][0]
        trp_y = trp_cd2_ce2_midpoints[trp_id][1]
        trp_z = trp_cd2_ce2_midpoints[trp_id][2]
        # iterate through each sugar residue
        # and each CH atom in that sugar
        for sugar_id, sugar_midpt in sugar_o5_c3_midpoints.iteritems():
            # pull out the sugar midpoint x, y, z
            sugar_x = sugar_o5_c3_midpoints[sugar_id][0]
            sugar_y = sugar_o5_c3_midpoints[sugar_id][1]
            sugar_z = sugar_o5_c3_midpoints[sugar_id][2]
            # calculate the distance to the TRP midpoint
            # sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
            dist = sqrt( (sugar_x - trp_x)**2 +
                         (sugar_y - trp_y)**2 +
                         (sugar_z - trp_z)**2 )
            # add the midpoint to CH distance to the list
            # for this current TRP's midpoint
            trp_midpt_to_sugar_midpt_dists[sugar_id] = dist
            # add the sugar midpoint xyz coords to the dict
            all_sugar_midpt_xyz_dict[sugar_id] = [sugar_x,
                                                  sugar_y,
                                                  sugar_z]

        # done checking all sugar midpoints for this specific TRP midpoint
        # determine which distance is the minimum
        # and update the data dictionary accordingly
        min_dist = None
        min_sugar_id = None
        for sugar_id, dist in trp_midpt_to_sugar_midpt_dists.iteritems():
            # if this is the first distance checked
            # or if this distance is smaller than it if it's not None
            # update the holders
            if (min_dist is None) or (dist < min_dist):
                min_dist = dist
                min_sugar_id = sugar_id

        # add this data to the dictionary
        # key: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain, value: distance
        # example: 3ACH_127_A+3ACH_Glc_310_A : 5.2
        trp_to_sugar_id = trp_id + "+" + min_sugar_id
        all_trp_midpt_to_sugar_midpt_dict[trp_to_sugar_id] = min_dist

    # we collected all the minimum TRP midpoints to sugar midpoints
    # ie) each TRP in the protein must have some minimum distance
    # to some sugar midpoint in the protein. that distance could be
    # large though. So go through the dictionary and prune out
    # TRP to sugar midpoint distances below some cutoff distance
    # (iterate through the same dictionary key and value pairs
    # and keep only values that are below the cutoff)
    # *** these are the TRP to sugar midpoint distances we want to futher analyze
    trp_midpt_to_sugar_midpt_dict = { trp_to_sugar_id:
                                          round(min_dist, 3)
                                      for trp_to_sugar_id, min_dist 
                                      in all_trp_midpt_to_sugar_midpt_dict.iteritems()
                                      if min_dist < trp_to_sugar_midpt_cutoff }
    # with the pruned TRP midpoint data, we can figure out which sugar midpoint was closest
    # and prune a dictionary of sugar midpoints to their xyz coords for later use
    # trp_to_sugar_id looks like: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain
    # so split on the "+" and pull out the sugar id to get its xyz coords
    sugar_midpt_xyz_dict = { trp_to_sugar_id.split("+")[-1]:
                                 all_sugar_midpt_xyz_dict[trp_to_sugar_id.split("+")[-1]]
                             for trp_to_sugar_id
                             in trp_midpt_to_sugar_midpt_dict.keys() }


###########################################
# CALCULATE PLANE BETWEEN TRP NE1 CD2 CH2 #
###########################################
    # arbitrary choice, but assuming this captures the
    # aromatic plane of the TRP residue
    # key: PDBname_TRPresnum_chain, value: [x, y, z]
    trp_ne1_cd2_ch2_planes = {}
    # trp_id = PDBname_TRPresnum_chain
    for trp_id, trp_res in trp_dict.iteritems():
        ne1_line = ""
        cd2_line = ""
        ch2_line = ""
        # get the NE1, CD2, and CH2 atom info
        for atom_line in trp_res:
            if atom_line[12:16].strip() == "NE1":
                # xyz info for NE1 atom
                ne1_line = atom_line
                ne1_x = float(atom_line[30:38].strip())
                ne1_y = float(atom_line[38:46].strip())
                ne1_z = float(atom_line[46:54].strip())
            if atom_line[12:16].strip() == "CD2":
                # xyz info for CD2 atom
                cd2_line = atom_line
                cd2_x = float(atom_line[30:38].strip())
                cd2_y = float(atom_line[38:46].strip())
                cd2_z = float(atom_line[46:54].strip())
            if atom_line[12:16].strip() == "CH2":
                # xyz info for CH2 atom
                ch2_line = atom_line
                ch2_x = float(atom_line[30:38].strip())
                ch2_y = float(atom_line[38:46].strip())
                ch2_z = float(atom_line[46:54].strip())
        # calculate the plane between the NE1, CD2, and CE2 atoms
        plane = Plane( Point3D(ne1_x, ne1_y, ne1_z),
                       Point3D(cd2_x, cd2_y, cd2_z),
                       Point3D(ch2_x, ch2_y, ch2_z) )
        # add to the dictionary
        trp_ne1_cd2_ch2_planes[trp_id] = plane


##########################################
# CALCULATE PLANE BETWEEN SUGAR C4 O5 C2 #
##########################################
    # arbitrary choice, but assuming this captures the
    # plane of the sugar residue
    # key: PDBname_SUGARresname_resnum_chain, value: [x, y, z]
    sugar_c4_o5_c2_planes = {}
    # sugar_id = PDBname_SUGARresname_resnum_chain
    for sugar_id, sugar_res in sugar_dict.iteritems():
        c4_line = ""
        o5_line = ""
        c2_line = ""
        # get the C4, O5, and C2 atom info
        for atom_line in sugar_res:
            if atom_line[12:16].strip() == "C4":
                # xyz info for C4 atom
                c4_line = atom_line
                c4_x = float(atom_line[30:38].strip())
                c4_y = float(atom_line[38:46].strip())
                c4_z = float(atom_line[46:54].strip())
            if atom_line[12:16].strip() == "O5":
                # xyz info for O5 atom
                o5_line = atom_line
                o5_x = float(atom_line[30:38].strip())
                o5_y = float(atom_line[38:46].strip())
                o5_z = float(atom_line[46:54].strip())
            if atom_line[12:16].strip() == "C2":
                # xyz info for C2 atom
                c2_line = atom_line
                c2_x = float(atom_line[30:38].strip())
                c2_y = float(atom_line[38:46].strip())
                c2_z = float(atom_line[46:54].strip())
        # calculate the plane between the C4, O5, and CE2 atoms
        plane = Plane( Point3D(c4_x, c4_y, c4_z),
                       Point3D(o5_x, o5_y, o5_z),
                       Point3D(c2_x, c2_y, c2_z) )
        # add to the dictionary
        sugar_c4_o5_c2_planes[sugar_id] = plane


################################################
# CALCULATE ANGLE BETWEEN TRP AND SUGAR PLANES #
################################################
    # the two planes calculated for the TRP and sugar residues
    # should roughly capture the CH-pi stacking between the residues
    # for each TRP residue that has a relevant calculation to make
    # i.e. it has a TRP midpoint to sugar midpoint distance below the cutoff
    # calculate the plane between the TRP and sugar residues
    # key: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain
    # value: TRP NE1, CD2, CH2 plane to sugar C4, O5, C2 plane angle
    # see below for why, but keeping the raw value and the value
    # you get after you subtract 180 and take the absolute value
    trp_to_sugar_plane_angle_raw_dict = {}
    trp_to_sugar_plane_angle_dict = {}
    for trp_to_sugar_id in trp_midpt_to_sugar_midpt_dict.keys():
        # pull out the trp_id
        trp_id = trp_to_sugar_id.split("+")[0]
        # pull out the sugar id too
        sugar_id = trp_to_sugar_id.split("+")[-1]

        # grab the Plane objects for the TRP and sugar residue
        trp_plane = trp_ne1_cd2_ch2_planes[trp_id]
        sugar_plane = sugar_c4_o5_c2_planes[sugar_id]

        # calculate the angle between the two planes
        # this function returns acos(some value) (this is radians)
        # so I need to finish the math
        # to get the angle in degrees, do: (180 * float(plane_angle)) / pi
        plane_angle_rads = trp_plane.angle_between(sugar_plane)
        plane_angle = round( (180 * float(plane_angle_rads)) / pi, 3)

        # add this raw plane_angle to the dictionary using the combined
        # trp and sugar id
        trp_to_sugar_plane_angle_raw_dict[trp_to_sugar_id] = plane_angle

        # ***** ASSUMPTION MADE HERE *****
        # we are going to assume that an angle of 0 is the equivalent
        # of an angle of 180. That is basically the sugar perfectly
        # stacked over the Trp in one conformation or flipped the other way
        # so we want our data to be adjusted so that is centered on one value
        # instead of being bimodal
        # we will therefore flip angle values that are 90 degrees or greater
        # by subtracting 180 from them and taking the absolute value
        # (angles are either 0 or some positive value)
        # basically, given the way that the planes are defined,
        # the difference is if the nitrogen of the TRP and the oxygen
        # of the sugar are on the same side or not
        # AND if the oxygen (and the two carbons) of the sugar
        # are pointing up or down with respect to the TRP
        if plane_angle >= 90:
            plane_angle = abs(plane_angle - 180)

        # add this plane angle to the dictionary using the combined
        # trp and sugar id
        # this could be the same as plane_angle_raw
        trp_to_sugar_plane_angle_dict[trp_to_sugar_id] = plane_angle
        


# TODO
# check that the distance between sugar O5 and C3 atoms is at least 2.8-2.9?
# calc that for all PDBs and see what the min and max is and if they are right
    '''
######################################################
# CALCULATE DIHEDRAL BETWEEN TRP CE2 CD2 SUGAR O5 C3 #
######################################################
    # for each TRP residue that has a relevant dihedral to calculate
    # i.e. it has a TRP midpoint to sugar midpoint distance below the cutoff
    # pull out the TRP residue from the trp_dict, get the relevant atoms
    # and then pull out the relevant atoms from the sugar_dict
    # key: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain, value: CE2-CD2-O5-C3 dihedral
    trp_to_sugar_dihedral_dict = {}
    for trp_to_sugar_id in trp_midpt_to_sugar_midpt_dict.keys():
        # pull out the trp_id
        trp_id = trp_to_sugar_id.split("+")[0]
        # pull out the sugar midpoint id too
        sugar_id = trp_to_sugar_id.split("+")[-1]

        # grab the TRP residue lines given the trp_id
        trp_res_lines = trp_dict[trp_id]
        # grab the sugar residue lines given the sugar_id
        sugar_res_lines = sugar_dict[sugar_id]

        # pull out the CE2 and CD2 atom coordinates for the TRP
        # start them out as None in case for some reason the TRP
        # residue doesn't have these atoms and we can't calculate the dihedral
        trp_CE2_xyz = None
        trp_CD2_xyz = None
        for trp_line in trp_res_lines:
             if trp_line[12:16].strip() == "CE2":
                trp_CE2_xyz = [float(trp_line[30:38].strip()),
                               float(trp_line[38:46].strip()),
                               float(trp_line[46:54].strip())]
             if trp_line[12:16].strip() == "CD2":
                 trp_CD2_xyz = [float(trp_line[30:38].strip()),
                                float(trp_line[38:46].strip()),
                                float(trp_line[46:54].strip())]

        # pull out the O5 and C3 atom coordinates for the sugar
        # start them out as None in case for some reason the sugar
        # residue doesn't have these atoms and we can't calculate the dihedral
        sugar_O5_xyz = None
        sugar_C3_xyz = None
        for sugar_line in sugar_res_lines:
             if sugar_line[12:16].strip() == "O5":
                sugar_O5_xyz = [float(sugar_line[30:38].strip()),
                               float(sugar_line[38:46].strip()),
                               float(sugar_line[46:54].strip())]
             if sugar_line[12:16].strip() == "C3":
                 sugar_C3_xyz = [float(sugar_line[30:38].strip()),
                                 float(sugar_line[38:46].strip()),
                                 float(sugar_line[46:54].strip())]

        # calculate the dihedral between CE2 - CD2 - O5 - C3
        trp_to_sugar_dihedral = round(calc_dihedral_angle(trp_CE2_xyz,
                                                          trp_CD2_xyz,
                                                          sugar_O5_xyz,
                                                          sugar_C3_xyz), 3)
        # add to the dictionary given the TRP to sugar id key
        trp_to_sugar_dihedral_dict[trp_to_sugar_id] = trp_to_sugar_dihedral
    '''


###############################################
# CALCULATE DIHEDRAL BETWEEN TRP CA CB CG CD1 #
###############################################
    # for all TRP residues that are CH-pi interacting based on the above criteria
    # calculate its chi2 value from atoms CA CB CG and CD1
    trp_chi2_dihedral_dict = {}
    for trp_to_sugar_id in trp_midpt_to_sugar_midpt_dict.keys():
        # pull out the TRPresnum_chain id
        trp_id = trp_to_sugar_id.split("+")[0]

        # grab the TRP residue lines given the trp_id
        trp_res_lines = trp_dict[trp_id]
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
        trp_chi2_dihedral_dict[trp_id] = trp_chi2_dihedral


###################
# DATA COLLECTION #
###################
    # TRP residue number and chain
    # pull from dictionary that has TRP and sugar id together
    # to ensure the order of collecting TRP and sugar ids corresponds
    # to the order that the TRP and sugar residues are actually interacting
    # key: PDBname_TRPresnum_chain+PDBname_SUGARresname_resnum_chain
    trp_to_sugar_ids = trp_midpt_to_sugar_midpt_dict.keys()
    trp_ids = [trp_to_sugar_id.split("+")[0]
               for trp_to_sugar_id in trp_to_sugar_ids]
    pdb_ids = [trp_id.split("_")[0]
                   for trp_id in trp_ids]
    # a pdb_id is something like 1abc
    # a pdb_name is something like 1abc_rosetta_suffix
    # a pdb_id was used to differentiate different data points collected
    # but, due to parsing (the _ specifically), the pdb_name wasn't used
    # but the pdb_name is desired because that's the actual filename
    # hence doing this little hack
    pdb_names = [pdb_name
                 for ii in range(len(pdb_ids))]
    trp_resnums = [int(trp_id.split("_")[1])
                   for trp_id in trp_ids]
    trp_chains = [trp_id.split("_")[2]
                  for trp_id in trp_ids]
    # sugar residue number, chain
    # sugar_id: SUGARresnum_chain_atomname
    sugar_ids = [trp_to_sugar_id.split("+")[1]
                 for trp_to_sugar_id in trp_to_sugar_ids]
    sugar_resnames = [sugar_id.split("_")[1]
                     for sugar_id in sugar_ids]
    sugar_resnums = [int(sugar_id.split("_")[2])
                     for sugar_id in sugar_ids]
    sugar_chains = [sugar_id.split("_")[3]
                    for sugar_id in sugar_ids]
    # TRP midpoint [x, y, z]
    trp_midpt_xyzs = [trp_cd2_ce2_midpoints[trp_id]
                      for trp_id in trp_ids]
    # sugar midpoint [x, y, z]
    sugar_midpt_xyzs = [sugar_o5_c3_midpoints[sugar_id]
                        for sugar_id in sugar_ids]
    # a little redundant, but explicit to ensure the correct order of things
    trp_midpt_to_sugar_midpt_dists = [float(trp_midpt_to_sugar_midpt_dict[trp_to_sugar_id])
                                      for trp_to_sugar_id in trp_to_sugar_ids]

    # TRP CE2 CD2 to sugar O5 C3 dihedrals
    # also a little redundant, but again, ensuring correct order
    #trp_to_sugar_dihedrals = [float(trp_to_sugar_dihedral_dict[trp_to_sugar_id])
    #                                for trp_to_sugar_id in trp_to_sugar_ids]

    # TRP NE1, CD2, CH2 plane to sugar C4, O5, C2 plane angle
    # also a little redundant, but again, ensuring correct order
    trp_to_sugar_plane_angles_raw = [float(trp_to_sugar_plane_angle_raw_dict[trp_to_sugar_id])
                                     for trp_to_sugar_id in trp_to_sugar_ids]
    trp_to_sugar_plane_angles = [float(trp_to_sugar_plane_angle_dict[trp_to_sugar_id])
                                 for trp_to_sugar_id in trp_to_sugar_ids]

    # TRP chi2 dihedrals
    trp_chi2_dihedrals = [float(trp_chi2_dihedral_dict[trp_id])
                          for trp_id in trp_ids]

    # all of the data lists at this point should be same length
    # so we'll just take one of those lists, pdb_ids
    # and use that as an iterator to then add each line of data
    # row by row to the dataframe
    # a dataframe is 0-indexed, so adding a row by doing len(df)
    # is good enough
    for ii in range(len(pdb_names)):
        df.loc[len(df)] = [pdb_names[ii],
                           trp_resnums[ii],
                           trp_chains[ii],
                           sugar_resnames[ii],
                           sugar_resnums[ii],
                           sugar_chains[ii],
                           trp_midpt_xyzs[ii],
                           sugar_midpt_xyzs[ii],
                           trp_midpt_to_sugar_midpt_dists[ii],
                           trp_to_sugar_plane_angles_raw[ii],
                           trp_to_sugar_plane_angles[ii],
                           #trp_to_sugar_dihedrals[ii],
                           trp_chi2_dihedrals[ii]]

# end time
end = time()
print "\nThis took %s seconds\n" %(round(end - start, 3))

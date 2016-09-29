#!/usr/bin/python
__author__="morganlnance"
__question__="ps3_q7"

'''
python disulfide_bond_detector.py <pdb filename>
'''

import sys, os


def calc_dist( x1, y1, z1, x2, y2, z2 ):
    from math import sqrt, pow
    return sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )


# grab commandline argument
try:
    pdb_filename = sys.argv[1]
except IndexError:
    print "\nPlease give me the path to a pdb file\n"
    sys.exit()

# check the file and open it
if not os.path.isfile( pdb_filename ):
    print "\nThe commandline argument you gave me does not seem to be a valid filepath\n"
    sys.exit()

try:
    with open( pdb_filename, "rb" ) as fh:
        pdb_file = fh.readlines()
except:
    print "\nI had an issue with opening your file ( %s )\n" %pdb_filename
    sys.exit()

# iterate through each ATOM line
sulfur_atom = "SG"
CYS_sulfur_lines = []
for line in pdb_file:
    # if ATOM line
    if line.startswith( "ATOM" ):
        # check if residue is a CYS and is not an alternate B conformation
        if line[17:20] == "CYS" and line[16:17] != 'B':
            # check if atom is a sulfur_atom
            if line[12:16].strip() == sulfur_atom:
                CYS_sulfur_lines.append( line.strip() )

# check the distance between all of the CYS sulfur atoms
# if less than 2.1A, then this is a disulfide
disulfide_bonds = []
for side1 in CYS_sulfur_lines:
    # get this side's unique name (chain_resnum)
    side1_chain = side1[21:22]
    side1_resnum = side1[22:26].strip()
    side1_uniq = '_'.join( [ side1_chain, side1_resnum ] )

    # get the xyz of this residue
    side1_x = float( side1[30:38].strip() )
    side1_y = float( side1[38:46].strip() )
    side1_z = float( side1[46:54].strip() )
    # check the xyz to every other residue
    for side2 in CYS_sulfur_lines:
        # get this side's unique name (chain_resnum)
        side2_chain = side2[21:22]
        side2_resnum = side2[22:26].strip()
        side2_uniq = '_'.join( [ side2_chain, side2_resnum ] )

        # if this isn't the same residue
        if side1_uniq != side2_uniq:
            side2_x = float( side2[30:38].strip() )
            side2_y = float( side2[38:46].strip() )
            side2_z = float( side2[46:54].strip() )

            # make sure the distance is less than 2.1
            if calc_dist( side1_x, side1_y, side1_z, side2_x, side2_y, side2_z ) <= 2.1:
                # make a unique name for this connection
                uniq_disulfide_way1 = '+'.join(  [ side1_uniq, side2_uniq ] )
                uniq_disulfide_way2 = '+'.join(  [ side2_uniq, side1_uniq ] )
                                
                # add to the disulfide connections list if not already there
                if uniq_disulfide_way1 not in disulfide_bonds and uniq_disulfide_way2 not in disulfide_bonds:
                    disulfide_bonds.append( uniq_disulfide_way1 )
                    disulfide_bonds.append( uniq_disulfide_way2 )
                    print "CYS %s %s  CYS %s %s" %( side1_chain, side1_resnum, side2_chain, side2_resnum )

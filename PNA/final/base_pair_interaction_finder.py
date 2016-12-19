#!/usr/bin/python
__author__="morganlnance"
__question__="final_q8"

'''
python base_pair_interaction_finder.py <pdb filename>
'''

import sys, os



G_atom_types = [ "O6", "N1", "N2" ]
C_atom_types = [ "O2", "N3", "N4" ]
A_atom_types = [ "N1", "N6" ]
U_atom_types = [ "O2", "O4", "N3" ]


def calc_dist( xyz_list1, xyz_list2 ):
    # pull the individual coordinates out from the list
    x1 = xyz_list1[ 0 ]
    y1 = xyz_list1[ 1 ]
    z1 = xyz_list1[ 2 ]
    x2 = xyz_list2[ 0 ]
    y2 = xyz_list2[ 1 ]
    z2 = xyz_list2[ 2 ]
    # calculate the distance
    from math import sqrt, pow
    return sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )

def get_xyz_coords( pdb_line ):
    return [ float( pdb_line[ 30:38 ].strip() ), 
             float( pdb_line[ 38:46 ].strip() ), 
             float( pdb_line[ 46:54 ].strip() ) ]

def get_res_name( pdb_line ):
    return pdb_line[17:20].strip()

def get_res_num( pdb_line ):
    return int( pdb_line[22:26] )

def get_res_chain( pdb_line ):
    return pdb_line[21:22]


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


# first get all of the ATOM lines
ATOM_lines = []
[ ATOM_lines.append( line ) for line in pdb_file if line.startswith( "ATOM" ) ]

# make a dictionary of residues
residue_dict = {}
for atom_line in ATOM_lines:
    uniq_name = '_'.join( [ str( get_res_num( atom_line ) ),
                            get_res_name( atom_line ),
                            get_res_chain( atom_line ) ] )
    if uniq_name not in residue_dict.keys():
        residue_dict[ uniq_name ] = []
    residue_dict[ uniq_name ].append( atom_line.strip() )

uniq_interactions = []
for uniq_name1 in residue_dict:
    for uniq_name2 in residue_dict:
        # if this is not the same residue
        if uniq_name1 != uniq_name2:
            # if these residues are not within +/- 1 residues of each other
            if int( uniq_name1.split( '_' )[0] ) != int( uniq_name2.split( '_' )[0] ) + 1:
                if int( uniq_name1.split( '_' )[0] ) != int( uniq_name2.split( '_' )[0] ) - 1:
                    # check only particular elements depending on residue type
                    atom1_lines = residue_dict[ uniq_name1 ]
                    atom2_lines = residue_dict[ uniq_name2 ]
                    for atom1 in atom1_lines:
                        res_type1 = atom1[17:20].strip()
                        atom_type1 = atom1[12:16].strip()
                        res_num1 = atom1[22:26].strip()
                        if ( res_type1 == 'G' and atom_type1 in G_atom_types ) \
                                or ( res_type1 == 'C' and atom_type1 in C_atom_types ) \
                                or ( res_type1 == 'A' and atom_type1 in A_atom_types ) \
                                or ( res_type1 == 'U' and atom_type1 in U_atom_types ):
                            for atom2 in atom2_lines:
                                res_type2 = atom2[17:20].strip()
                                atom_type2 = atom2[12:16].strip()
                                res_num2 = atom2[22:26].strip()
                                if ( res_type2 == 'G' and atom_type2 in G_atom_types ) \
                                        or ( res_type2 == 'C' and atom_type2 in C_atom_types ) \
                                        or ( res_type2 == 'A' and atom_type2 in A_atom_types ) \
                                        or ( res_type2 == 'U' and atom_type2 in U_atom_types ):
                                    if calc_dist( get_xyz_coords( atom1 ), get_xyz_coords( atom2 ) ) <= 3.0:
                                        # for use in checking in pymol
                                        #print "dist /1z43//A/%s`%s/%s, /1z43//A/%s`%s/%s" %( res_type1, res_num1, atom_type1, res_type2, res_num2, atom_type2 )

                                        uniq_interaction = uniq_name1 + '+' + uniq_name2
                                        rev_uniq_interaction = uniq_name2 + '+' + uniq_name1
                                        if uniq_interaction not in uniq_interactions:
                                            if rev_uniq_interaction not in uniq_interactions:
                                                uniq_interactions.append( uniq_interaction )
                                                uniq_interactions.append( rev_uniq_interaction )
                                                break

print "\nBase pair interactions for %s\n" %pdb_filename.split( '/' )[-1]
uniq_interactions.sort()
for interaction in uniq_interactions:
    interaction1 = interaction.split( '+' )[0]
    interaction2 = interaction.split( '+' )[1]
    res_type1 = interaction1.split( '_' )[1]
    res_num1 = interaction1.split( '_' )[0]
    res_type2 = interaction2.split( '_' )[1]
    res_num2 = interaction2.split( '_' )[0]
    if res_type1 == 'G' and res_type2 == 'C':
        print "%s %s H-bonds to %s %s Watson-Crick" %( res_type1, res_num1, res_type2, res_num2 )
    elif res_type1 == 'C' and res_type2 == 'G':
        print "%s %s H-bonds to %s %s Watson-Crick" %( res_type1, res_num1, res_type2, res_num2 )
    elif res_type1 == 'A' and res_type2 == 'U':
        print "%s %s H-bonds to %s %s Watson-Crick" %( res_type1, res_num1, res_type2, res_num2 )
    elif res_type1 == 'U' and res_type2 == 'A':
        print "%s %s H-bonds to %s %s Watson-Crick" %( res_type1, res_num1, res_type2, res_num2 )
    else:
        print "%s %s H-bonds to %s %s non Watson-Crick" %( res_type1, res_num1, res_type2, res_num2 )
        
print

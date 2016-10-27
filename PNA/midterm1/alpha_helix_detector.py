#!/usr/bin/python
__author__="morganlnance"
__question__="midterm1_8"

'''
python alpha_helix_detector.py <pdb filename>
'''

import sys, os


# for use in converting to FASTA format
AA_name3_to_name1 = { "ALA":'A', 
                      "CYS":'C', 
                      "ASP":'D', 
                      "GLU":'E', 
                      "PHE":'F', 
                      "GLY":'G', 
                      "HIS":'H', 
                      "ILE":'I', 
                      "LYS":'K', 
                      "LEU":'L', 
                      "MET":'M', 
                      "ASN":'N', 
                      "PRO":'P', 
                      "GLN":'Q', 
                      "ARG":'R', 
                      "SER":'S', 
                      "THR":'T', 
                      "VAL":'V', 
                      "TRP":'W', 
                      "TYR":'Y' }


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

def get_atom_name( pdb_line ):
    return pdb_line[ 12:16 ].strip()


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

# then break up the atom lines into a dictionary based on residue number
residue_lines = {}
num_to_uniq_res_name = {}
ii = -1  # first residue in pdb file will get number 0
for line in ATOM_lines:
    try:
        uniq_res_name = '_'.join( [ AA_name3_to_name1[ line[17:20] ], # single-letter residue name
                                    line[21:22], # residue chain
                                    line[22:26].strip() ] ) # residue number
    # skip non-canonical amino acids
    except KeyError:
        continue

    # check the current residue to the last unique residue
    # if this is a new residue, create a new entry in the dict
    if uniq_res_name not in residue_lines.keys():
        residue_lines[ uniq_res_name ] = [ line ]
        ii += 1
        num_to_uniq_res_name[ ii ] = uniq_res_name
    # else, we've seen this residue just last time, so add the next atom to it
    else:
        residue_lines[ uniq_res_name ].append( line )

# get the protein length from the number of unique residues
protein_len = len( residue_lines.keys() )

# for each unique residue
H_res_nums = []  # residues to be designated as a helix
for ii in range( protein_len ):
    # if this isn't the last three residues
    if ii + 4 < protein_len:
        # get the current residue and the ii + 4 residue
        res1_lines = residue_lines[ num_to_uniq_res_name[ ii ] ]
        res2_lines = residue_lines[ num_to_uniq_res_name[ ii + 4 ] ]

        # grab the xyz coordinates of the backbone O of residue ii
        # backbone O is only the letter O
        for line in res1_lines:
            if get_atom_name( line ) == 'O':
                backbone_O_xyz = get_xyz_coords( line )
                break

        # grab the xyz coordinates of the backbone N of residue ii + 4
        # backbone N is only the letter N
        for line in res2_lines:
            if get_atom_name( line ) == 'N':
                backbone_N_xyz = get_xyz_coords( line )
                break

        # check the distance between all of the backbone N and O atoms
        # if less than 4A, then this is a hydrogen bond
        if calc_dist( backbone_O_xyz, backbone_N_xyz ) <= 4.0:
            # residue ii and ii + 4 should have an 'H' designation
            if ii not in H_res_nums:
                H_res_nums.append( ii )
            if ii + 4 not in H_res_nums:
                H_res_nums.append( ii + 4 )

# for each unique residue name, pull out the single-letter amino acid code and give it an 'H' or '-' designation
# holds the amino acide single-letter code for each residue
# my unique name was <res name>_<res chain>_<res num>
FASTA_aa = [ num_to_uniq_res_name[ ii ].split( '_' )[0] for ii in range( protein_len ) ]
# holds the secondary structure code for each residue
FASTA_ss = []
for ii in range( protein_len ):
    if ii in H_res_nums:
        FASTA_ss.append( 'H' )
    else:
        FASTA_ss.append( '-' )

# print the FASTA secondary structure designation above the corresponding residue
# write the same data to a file
pdb_name = pdb_filename.split( '/' )[-1].split( ".pdb" )[0]
FASTA_len = 80
print pdb_name, "FASTA sequence describing alpha helical secondary structure\n"
with open( pdb_name + ".fasta", "wb" ) as fh:
    for ii in range( 0, protein_len, FASTA_len ):
        if protein_len - ii < FASTA_len:
            print ''.join( FASTA_ss[ ii : protein_len ] )
            print ''.join( FASTA_aa[ ii : protein_len ] )
            # write to file
            fh.write( ''.join( FASTA_ss[ ii : protein_len ] ) )
            fh.write( "\n" )
            fh.write( ''.join( FASTA_aa[ ii : protein_len ] ) )
            fh.write( "\n" )
        else:
            print ''.join( FASTA_ss[ ii : ii + FASTA_len ] )
            print ''.join( FASTA_aa[ ii : ii + FASTA_len ] )
            # write to file
            fh.write( ''.join( FASTA_ss[ ii : ii + FASTA_len ] ) )
            fh.write( "\n" )
            fh.write( ''.join( FASTA_aa[ ii : ii + FASTA_len ] ) )
            fh.write( "\n" )

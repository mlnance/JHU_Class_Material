#!/usr/bin/python
__author__="morganlnance"


import numpy as np
from math import sin, cos, pi
import sys


def to_radians( degrees ):
    # degrees * ( pi / 180 degrees ) = radians
    return degrees * ( pi / 180.0 )


# open the input .mop file
try:
    in_file = sys.argv[1]
except IndexError:
    print "\nYou did not give me an input file. I need a .mop file.\n"
    sys.exit()
try:
    with open( in_file, "r" ) as fh:
        lines = fh.readlines()
except:
    print "\nThere was something wrong with your input file.\n"
    sys.exit()

# insantiate an iterator so we can count atoms
atom_num = 1
res_num = 1
# to hold lines converted to PDB format
pdb_lines = []
occupancy = "%.2f" %1.00
t_factor = "%.2f" %30.00

# go through line by line to pull out the necessary information
# we need bond length, theta, and phi to calculate Cartesian x, y, z
for line in lines:
    # split the columns
    split_line = line.split()
    # .mop file column information
    # atom_name, b_length, 1, theta, 1, phi, 1, con1, con2, con3
    atom_name = split_line[0]
    b_length = float( split_line[1] )
    theta = float( split_line[3] )
    phi = float( split_line[5] )

    # get the values needed to compute the x, y, z matrix
    # see PowerPoint for information
    # matrix: [ [ 00, 01, 02 ], [ 10, 11, 12 ], [ 20, 21, 22 ] ]
    pos00 = -( cos( to_radians( theta ) ) )
    pos01 = ( sin( to_radians( theta ) ) )
    pos02 = 0.0
    pos10 = -( cos( to_radians( phi ) ) ) * ( sin( to_radians( theta) ) )
    pos11 = -( cos( to_radians( phi ) ) ) * -( cos( to_radians( theta ) ) )
    pos12 = -( sin( to_radians( phi) ) )
    pos20 = -( sin( to_radians( phi ) ) ) * ( sin( to_radians( theta) ) )
    pos21 = -( sin( to_radians( phi ) ) ) * ( cos( to_radians( theta ) ) )
    pos22 = ( cos ( to_radians( phi ) ) )

    # set up the rotation matrix
    rotation_matrix = np.matrix( [ [pos00, pos01, pos02 ], 
                                   [pos10, pos11, pos12 ], 
                                   [pos20, pos21, pos22 ] ] )
    # set up the coordinate vector
    coordinate_vector = np.matrix( [ [ b_length ], [ 0 ], [ 0 ] ] )

    # mutliply the two matrices to get the x, y, z coordinates
    xyz_coords = rotation_matrix * coordinate_vector
    x = "%.3f" %xyz_coords[0][0]
    y = "%.3f" %xyz_coords[1][0]
    z = "%.3f" %xyz_coords[2][0]

    ## prepare to write .pdb line
    pdb_line = "ATOM  " + \
        ' ' * ( 5 - len( str( atom_num ) ) ) + str( atom_num ) + \
        ' ' * ( 4 - len( atom_name ) ) + atom_name + \
        ' ' + \
        " INB " + \
        'A' + \
        ' ' * ( 4 - len( str( res_num ) ) ) + str( res_num ) + \
        ' ' * 4 + \
        ' ' * ( 8 - len( x ) ) + x + \
        ' ' * ( 8 - len( y ) ) + y + \
        ' ' * ( 8 - len( z ) ) + z + \
        ' ' * ( 6 - len( occupancy ) ) + occupancy + \
        ' ' * ( 6 - len( t_factor ) ) + t_factor + \
        ' ' * 11 + \
        atom_name + \
        '   \n'
    pdb_lines.append( pdb_line )
    
    # atom number increase
    atom_num += 1


# write the PDB lines to a file
with open( "out.pdb", "w" ) as fh:
    fh.writelines( pdb_lines )

#!/usr/bin/python
__author__="morganlnance"
__question__="ps11_q7"

'''
python 1t7p_interaction_script.py <pdb filename>
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


# class object for holding residue details of a particular atom
class AtomLine:
    def __init__( self, line ):
        self.line = line.strip()
        self.atom_num = int( self.line[6:11].replace( ' ', '' ) )
        self.atom_name = str( self.line[12:16].replace( ' ', '' ) )
        self.res_name = str( self.line[17:20].replace( ' ', '' ) )
        self.res_chain = str( self.line[21:22] )
        self.res_num = int( self.line[22:26].replace( ' ', '' ) )
        self.x = float( self.line[30:38].replace( ' ', '' ) )
        self.y = float( self.line[38:46].replace( ' ', '' ) )
        self.z = float( self.line[46:54].replace( ' ', '' ) )
        self.element = str( self.line[76:78].replace( ' ', '' ) )
        self.charge = str( self.line[78:80].replace( ' ', '' ) )
        

# first, pull out all of the lines and coordinates for all heavy atoms of DG3
DG3_lines = []
for line in pdb_file:
    if line.startswith( "HETATM" ):
        atom_line = AtomLine( line.strip() )
        if atom_line.res_name == "DG3":
            # skip hydrogens
            if not atom_line.element == 'H':
                DG3_lines.append( atom_line )
# DEBUG to ensure you just pulled out DG3
#print "\n".join( [ line.line for line in DG3_lines ] )

# now find all C to C and polar atom to polar atom contacts between
# DG3 and the rest of the structure
polar_atoms = [ 'N', 'O', 'P', 'MG' ]  # assuming these can H-bond or do something interesting
pymol_dist_commands = []
vdw_counter = 1
polar_counter = 1
for line in pdb_file:
    if line.startswith( "ATOM" ) or line.startswith( "HETATM" ):
        atom_line = AtomLine( line.strip() )
        # skip the DG3 residue because we don't want self-self contacts
        if not atom_line.res_name == "DG3":
            # skip hydrogens
            if not atom_line.element == 'H':
                # skip waters
                if not atom_line.res_name == "HOH":
                    # check distance of this elemnt to each element in DG3
                    # wow this will take forever
                    for DG3 in DG3_lines:
                        # dist ___, /1t7p//chain/DG3`resnum/atomname, /1t7p//chain/resname`resnum/atomname
                        if calc_dist( DG3.x, DG3.y, DG3.z, atom_line.x, atom_line.y, atom_line.z ) <= 3.5:
                            # if carbon to carbon, mark blue as Van der Waals
                            if DG3.element == 'C' and atom_line.element == 'C':
                                pymol_dist_commands.append( "dist vdw%s, /1t7p//%s/%s`%s/%s, /1t7p//%s/%s`%s/%s\n" %( str( vdw_counter ), 
                                                                                                                      DG3.res_chain, 
                                                                                                                      DG3.res_name, 
                                                                                                                      DG3.res_num, 
                                                                                                                      DG3.atom_name, 
                                                                                                                      atom_line.res_chain, 
                                                                                                                      atom_line.res_name, 
                                                                                                                      atom_line.res_num, 
                                                                                                                      atom_line.atom_name ) )
                                vdw_counter += 1
                            # if some kind of polar interaction
                            elif DG3.element in polar_atoms and atom_line.element in polar_atoms:
                                pymol_dist_commands.append( "dist polar%s, /1t7p//%s/%s`%s/%s, /1t7p//%s/%s`%s/%s\n" %( str( polar_counter ), 
                                                                                                                        DG3.res_chain, 
                                                                                                                        DG3.res_name, 
                                                                                                                        DG3.res_num, 
                                                                                                                        DG3.atom_name, 
                                                                                                                        atom_line.res_chain, 
                                                                                                                        atom_line.res_name, 
                                                                                                                        atom_line.res_num, 
                                                                                                                        atom_line.atom_name ) )
                                polar_counter += 1

# write commands to a pymol script
cwd = os.getcwd() + '/'
pml_filename = cwd + "1t7p_interactions.pml"
with open( pml_filename, "wb" ) as fh:
    fh.write( "load 1t7p.pdb\n" )
    fh.write( "hide everything, resn hoh\n" )
    fh.writelines( pymol_dist_commands )
    fh.write( "group vdw_contacts, vdw*\n" )
    fh.write( "group polar_contacts, polar*\n" )
    fh.write( "show sticks, resn DG3\n" )
    fh.write( "color magenta, vdw_contacts\n" )
    fh.write( "color yellow, polar_contacts\n" )
    fh.write( "hide labels\n" )
    fh.write( "show spheres, resn MG\n" )
    fh.write( "set sphere_scale, 0.4\n" )
    fh.write( "bg_color white\n" )
    fh.write( "set_view (\
     0.059504092,    0.576347291,    0.815033674,\
    -0.203007400,    0.806406200,   -0.555421352,\
    -0.977367759,   -0.132408097,    0.164987102,\
     0.000000000,    0.000000000,  -49.289268494,\
    45.421001434,   24.983999252,    0.598999977,\
    41.054149628,   61.524391174,  -20.000000000 )" )
print "\nDumped '1t7p_interactions.pml' in your working directory. Run with @1t7p_interactions.pml in pymol"
print "Van der Waals contacts are colored in magenta"
print "H-bond/electrostatic contacts are colored in yellow\n"

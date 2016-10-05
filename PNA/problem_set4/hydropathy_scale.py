#!/usr/bin/python
__author__="morganlnance"
__question__="ps4_q7"


'''
python hydropathy_scale.py <pdb file>
'''

# imports
import sys, os
import matplotlib
matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt


# take in the argument from the commandline (should be a pdb file)
try:
    pdb_file = sys.argv[1]
except:
    print "\nPlease give me a pdb file\n"
    sys.exit()
# check that the pdb file exists and ends with .pdb
if not os.path.isfile( pdb_file ):
    print "\nWhatever file you gave me is not actually a file. Check your input\n"
    sys.exit()
if not pdb_file.endswith( ".pdb" ):
    print "\nAre you sure you gave me a pdb file? Your file didn't end with '.pdb'\n"
    sys.exit()


# global variables
# adapted from the Wimley-White hydrophobicity scale
hydrophobicity_scale = { 'ALA' : 0.17, 
                         'ARG' : 0.81, 
                         'ASN' : 0.42, 
                         'ASP' : -0.07, 
                         'CYS' : -0.24, 
                         'GLN' : 0.58, 
                         'GLU' : -0.01, 
                         'GLY' : 0.01, 
                         'HIS' : 0.17, 
                         'ILE' : -0.31, 
                         'LEU' : -0.56, 
                         'LYS' : 0.99, 
                         'MET' : -0.23, 
                         'PHE' : -1.13, 
                         'PRO' : 0.45, 
                         'SER' : 0.13, 
                         'THR' : 0.14, 
                         'TRP' : -1.85, 
                         'TYR' : -0.94, 
                         'VAL' : 0.07 }
window_size = 19 # should be odd!
middle_of_window = ( window_size / 2 ) + 1  # ex. 3/2 in python is 1, so +1 gives 2 which is in the middle of 1 and 3


# open and read the pdb
try:
    with open( pdb_file, "rb" ) as fh:
        pdb_contents = fh.readlines()
except:
    print "\nI couldn't open the file you gave me for some reason :(\n"
    sys.exit()

# get the number of residues in the pdb
uniq_res_names = []
for line in pdb_contents:
    if line.startswith( "ATOM" ):
        # uniq_name is resname_reschain_resnum
        uniq_name = line[17:20] + '_' + line[21:22] + '_' + line[22:26].strip()
        if uniq_name not in uniq_res_names:
            uniq_res_names.append( uniq_name )

# cycle through the residue names and look at the middle of the residue in chunks of <window_size>
resnums = []
hydrophobicities = []
for ii in range( len( uniq_res_names ) ):
    # splitting on '_' because that's how I created my unique name
    if ii+window_size <= len( uniq_res_names ):
        # get the hydrophobicity of this residue in the middle of the window
        res_in_mid_of_window, resnum_in_mid_of_window = [ ( res.split( '_' )[0], res.split( '_' )[2] ) for res in uniq_res_names[ ii : ii+window_size ] ][ middle_of_window ]
        try:
            hydrophobicity = hydrophobicity_scale[ res_in_mid_of_window ]
            resnums.append( resnum_in_mid_of_window )
            hydrophobicities.append( hydrophobicity )
        except:
            # skip residues that are not canonical
            continue


# plot the figure
plt.figure(figsize=(16,12))
plt.plot( resnums, hydrophobicities )
plt.title( "Wimley-White Hydrophobicity of Residues Found in %s" %pdb_file )
plt.xlabel( "Residue Number" )
plt.ylabel( "Wimley-White Hydrophobicity" )
plt.show()
#plt.close()

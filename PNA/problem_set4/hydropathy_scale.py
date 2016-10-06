#!/usr/bin/python
__author__="morganlnance"
__question__="ps4_q7"

'''
> python hydropathy_scale.py <pdb file>
This will dump a .png file of the hydropathy plot in your current working directory
This script uses the octanol values from Wimley-White
'''

# imports
import sys, os
#import matplotlib
#matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 22 } )


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
# adapted from the Wimley-White octanol hydrophobicity scale
hydrophobicity_scale = { 'ALA' : 0.50, 
                         'ARG' : 1.81, # charged
                         'ASN' : 0.85, 
                         'ASP' : 3.64,  # charged
                         'CYS' : -0.02, 
                         'GLN' : 0.77, 
                         'GLU' : 3.63, # charged
                         'GLY' : 1.15, 
                         'HIS' : 0.11, 
                         'ILE' : -1.12, 
                         'LEU' : -1.25, 
                         'LYS' : 2.80, # charged
                         'MET' : -0.67, 
                         'PHE' : -1.71, 
                         'PRO' : 0.14, 
                         'SER' : 0.46, 
                         'THR' : 0.25, 
                         'TRP' : -2.09, 
                         'TYR' : -0.71, 
                         'VAL' : -0.46 }
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
        # get the hydrophobicity of the window (need resname from first part of uniq_name )
        resnames_in_window = [ res.split( '_' )[0] for res in uniq_res_names[ ii : ii+window_size ] ]
        # pull the window out again and get the residue number that is in the middle of the window (need resnum from last part of uniq_name
        resnum_in_mid_of_window = [ res.split( '_' )[-1] for res in uniq_res_names[ ii : ii+window_size ] ][ middle_of_window ]
        # get the average hydrophobicity for this window
        hydrophobicity = sum( [ hydrophobicity_scale[ resname ] for resname in resnames_in_window ] ) / len( resnames_in_window )
        # append the residue number and the hydrophobicity
        resnums.append( resnum_in_mid_of_window )
        hydrophobicities.append( hydrophobicity )


# plot the figure
plt.figure(figsize=(20,14))
plt.plot( resnums, hydrophobicities )
plot_title = "Wimley-White_Octanol_Hydrophobicity_of_Residues_in_%s" %pdb_file.split( ".pdb" )[0]
plt.title( plot_title )
plt.xlabel( "Residue Number (Window size of %s)" %window_size )
plt.ylabel( "Wimley-White Octanol Hydrophobicity" )
plt.savefig( plot_title, dpi=120, transparent=True )
#plt.show()
#plt.close()

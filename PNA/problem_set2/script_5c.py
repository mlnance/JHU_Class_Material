#!/usr/bin/python
__author__="morganlnance"
__question__="ps2_5c"


from pymol import cmd

'''
Question: use logic to select buried waters? Or do so explicitly?
'''

# load
cmd.load("1ij1.pdb")

# pretty-fy
cmd.show( "cartoon" )
cmd.hide( "lines" )

# show only the buried waters as spheres, otherwise remove the waters
# C204, B203, A202, A203
cmd.select( "delete_waters", "resn hoh and not ( resi 203+202 and chain A + resi 203 and chain B + resi 204 and chain C )" )
cmd.remove( "delete_waters" )
cmd.delete( "delete_waters" )
cmd.select( "buried_waters", "resn hoh" )
cmd.show( "spheres", "buried_waters" )
cmd.set( "sphere_scale", 0.5 )
cmd.h_add( "buried_waters" )
cmd.deselect()

# select the Thr residues around the buried waters
cmd.select( "thr_residues", "resi 12 and chain A + resi 12 and chain B + resi 12 and chain C" )
cmd.show( "sticks", "thr_residues and !(name c+o+n )" )
cmd.deselect()

# get the distances between the O of HOH and O of the 3 Thr residues
cmd.dist( "thr_hoh_A", "/1ij1//A/THR`12/OG1", "/1ij1//A/HOH`202/O" )
cmd.dist( "thr_hoh_B", "/1ij1//B/THR`12/OG1", "/1ij1//B/HOH`203/O" )
cmd.dist( "thr_hoh_C", "/1ij1//C/THR`12/OG1", "/1ij1//C/HOH`204/O" )

# set the view
cmd.set_view ( "0.808693528, 0.356420070, 0.467955768,\
                   -0.388590366, 0.920925915, -0.029888431,\
                   -0.441604972, -0.157672778, 0.883246481,\
                   0.000000000, 0.000000000, -50.762001038,\
                   -2.018598557, -0.170596123, 11.568523407,\
                   15.002751350, 86.521240234, -20.000000000" )

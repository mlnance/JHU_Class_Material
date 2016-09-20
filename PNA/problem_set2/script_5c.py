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

# show waters as sticks
cmd.show( "spheres", "resn hoh" )
cmd.set( "sphere_scale", 0.5 )
cmd.h_add( "resn hoh" )

# select buried waters
# C204, B203, A203, A202

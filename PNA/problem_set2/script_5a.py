#!/usr/bin/python
__author__="morganlnance"
__question__="ps2_5a"


from pymol import cmd

# load
cmd.load("1gcl.pdb")
cmd.load("2zta.pdb")

# pretty-fy
cmd.remove( "resn hoh" )
cmd.show( "cartoon" )
cmd.hide( "lines" )

# align on chain A
cmd.super( "1gcl and chain A", "2zta and chain A" )

# show CA atoms
cmd.show( "spheres", "name ca" )
cmd.set( "sphere_scale", "0.35" )


# set the final view
#cmd.set_view (\
#    0.298093528,    0.497441351,    0.814672232,\
#        -0.506814241,    0.805717707,   -0.306527108,\
#        -0.808874786,   -0.321513742,    0.492289960,\
#        -0.000013085,   -0.000000924,  -96.456314087,\
#        22.441783905,   24.351091385,    6.097001553,\
#        63.853878021,  129.056365967,  -20.000000000 )

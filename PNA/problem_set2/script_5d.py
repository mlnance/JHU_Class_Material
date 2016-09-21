#!/usr/bin/python
__author__="morganlnance"
__question__="ps2_5d"


from pymol import cmd, util


# to allow variables to hold additional information
class my_obj:
    pass

# load
cmd.load("216l.pdb")

# clean: remove hoh, color by chain
cmd.remove( "resn hoh" )
util.cbc()

# create the symmetry mates
cmd.symexp( "sym", "216l", "(216l)", 4 )
cmd.show( "cartoon" )
cmd.hide( "lines" )

# show the Trp at resi 44
cmd.select( "trp_44", "resi 44 and !( name c+o+n )" )
cmd.show( "sticks", "trp_44" )
cmd.color( "magenta", "trp_44" )
cmd.dist( "trp_trp_dist", "/sym010000-1//A/TRP`44/CZ3", "/216l//A/TRP`44/CZ3" )
cmd.deselect()

# set the view
cmd.set_view('''\
    0.104914360,   -0.982791662,    0.152021840,\
        -0.108071871,   -0.163226306,   -0.980651259,\
        0.988590658,    0.086455040,   -0.123336211,\
        -0.000002299,    0.000189751, -121.594604492,\
        56.129817963,   13.524446487,   -1.704995871,\
        73.446876526,  169.744995117,  -20.000000000''' )

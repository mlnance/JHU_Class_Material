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
cmd.deselect()

# set the view
cmd.set_view ('''\
    -0.850933135,   -0.424174100,    0.309816211,\
         -0.236527279,   -0.217213050,   -0.947033763,\
         0.469003528,   -0.879143298,    0.084505856,\
         0.000206053,    0.000086367, -139.347091675,\
         62.938968658,    9.275196075,    1.290278196,\
         94.235290527,  184.458694458,  -20.000000000''' )

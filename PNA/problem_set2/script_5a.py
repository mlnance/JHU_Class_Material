#!/usr/bin/python
__author__="morganlnance"
__question__="ps2_5a"
# run with pymol -r script_5a.py


from pymol import cmd


# to allow variables to hold additional information
class my_obj:
    pass

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

# select a-position residues on 2zta (dimer)
_2zta = my_obj()
_2zta.name = "2zta_a_residues"
_2zta.selection = "resi " + '+'.join( str(ii) for ii in range( 2, 31, 7 ) ) + " and 2zta"
cmd.select( _2zta.name, _2zta.selection )
cmd.show( "sticks", _2zta.name )
cmd.deselect()

# select a-position residues on 1gcl (dimer)
_1gcl = my_obj()
_1gcl.name = "1gcl_a_residues"
_1gcl.selection = "resi " + '+'.join( str(ii) for ii in range( 2, 31, 7 ) ) + " and 1gcl"
cmd.select( _1gcl.name, _1gcl.selection )
cmd.show( "sticks", _1gcl.name )
cmd.deselect()
    


# set the final view
cmd.set_view (\
    "0.133721262,    0.506625175,    0.851732850,\
        -0.829290450,    0.527751088,   -0.183717906,\
        -0.542578638,   -0.681767941,    0.490712464,\
        0.000017121,   -0.000036474,  -90.170852661,\
        19.357597351,   25.288087845,    5.837380409,\
        65.057907104,  115.283477783,  -20.000000000" )

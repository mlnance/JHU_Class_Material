#!/usr/bin/python
__author__="morganlnance"
__question__="ps2_5b"


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

# select d-position residues on 2zta (dimer)
_2zta = my_obj()
_2zta.name = "2zta_d_residues"
_2zta.selection = "resi " + '+'.join( str( ii ) for ii in range( 5, 31, 7 ) ) + " and 2zta"
cmd.select( _2zta.name, _2zta.selection )
cmd.show( "sticks", _2zta.name )
cmd.deselect()

# select d-position residues on 1gcl (dimer)
_1gcl = my_obj()
_1gcl.name = "1gcl_d_residues"
_1gcl.selection = "resi " + '+'.join( str( ii ) for ii in range( 5, 31, 7 ) ) + " and 1gcl"
cmd.select( _1gcl.name, _1gcl.selection )
cmd.show( "sticks", _1gcl.name )
cmd.deselect()
    


# set the final view
cmd.set_view ('''\
     0.054921620,    0.551247716,    0.832531750,\
    -0.880743742,    0.419553339,   -0.219698429,\
    -0.470399141,   -0.721180677,    0.508550465,\
    -0.000057466,    0.000006534,  -72.649322510,\
    19.602827072,   26.242700577,    7.284172058,\
    35.317283630,  109.982116699,  -20.000000000''' )

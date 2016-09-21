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
cmd.align( "1gcl and chain A and resi 9", "2zta and chain A and resi 9" )
#cmd.super( "1gcl and chain A", "2zta and chain A" )

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
cmd.set_view ('''\
    0.477873474,    0.178628594,    0.860074282,\
        -0.098605908,    0.983823717,   -0.149543911,\
        -0.872874737,   -0.013345651,    0.487757653,\
        0.000025588,   -0.000082150,  -67.047599792,\
        24.529411316,   18.095930099,   14.328298569,\
         62.294055939,   86.796539307,  -20.000000000''' )

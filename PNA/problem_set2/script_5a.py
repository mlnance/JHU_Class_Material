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
#cmd.show( "spheres", "name ca" )
#cmd.set( "sphere_scale", "0.35" )

# change cartoon transparency
cmd.set( "cartoon_transparency", "0.5" )

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
     0.532562733,    0.395187587,    0.748467803,\
    -0.627743840,    0.777577996,    0.036103889,\
    -0.567724824,   -0.489071995,    0.662189424,\
     0.000053290,   -0.000126574,  -48.477600098,\
    24.474441528,   18.028413773,   14.390577316,\
    43.106315613,   89.651206970,  -20.000000000''' )

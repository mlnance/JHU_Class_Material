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

# align structures
cmd.align( "1gcl and chain A", "2zta and chain A" )


# show CA atoms
#cmd.show( "spheres", "name ca" )
#cmd.set( "sphere_scale", "0.35" )

# set cartoon transparency
cmd.set( "cartoon_transparency", "0.5" )

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
     0.319376051,    0.492265671,    0.809733808,\
    -0.882474005,    0.465868622,    0.064847313,\
    -0.345308959,   -0.735281110,    0.583199799,\
    -0.000008503,   -0.000050932,  -71.901115417,\
    17.041526794,   16.113079071,    8.648363113,\
    34.564842224,  109.229652405,  -20.000000000''' )

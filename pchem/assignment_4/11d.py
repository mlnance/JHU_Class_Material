#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )
from math import exp


dH = 418  # kJ/mol
dS = 1.32 # kJ/( mol * K )
R = 0.008314 # kJ/( mol * K )


def dG( dH, T, dS ):
    return dH - ( T * dS )


def Keq( dH, T, dS ):
    return exp( -dG( dH, T, dS ) / ( R * T ) )

def fn( dH, T, dS ):
    return 1 / ( 1 + Keq( dH, T, dS ) )

# fn = 1 / 1 + Keq and fd = Keq / 1 + Keq
# Keq = Kd in the N <--> D reaction

T_data = []
fn_data = []
for T in range( 280, 341 ):
    T_data.append( T )
    fn_data.append( fn( dH, T, dS ) )

plt.plot( T_data, fn_data )
plt_title = "Generic Denaturation Curve"
plt.title( plt_title )
plt.xlabel( "T in Kelvin" )
plt.ylabel( "Fraction native" )
plt.savefig( plt_title.replace( ' ', '_' ), dpi=120, transparent=True )
plt.close()


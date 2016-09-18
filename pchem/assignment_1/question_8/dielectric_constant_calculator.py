#!/usr/bin/python
__author__="morganlnance"


import matplotlib
matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt
from math import log


def dielectric( T ):
    """
    Calculate the dielectric constant using the empirical formula and a temperature
    :param T: float( temperature )
    :return float
    """
    e = 87.740 - ( 0.40008*T ) + ( 0.0009398 * (T**2) ) - ( 0.000001410 * (T**3) )

    return e


# T 0 through 100 Celsius
x_val = range( 101 )

# T 0 through 100
y_val = [ dielectric( T ) for T in x_val ]
plt.plot( x_val, y_val )
plt.xlabel( "T in degrees C" )
plt.ylabel( "e (dielectric constant of water)" )
plt.title( "Temperature dependence of dielectric constant of water" )


# save plot
plot_title = "e (dielectric of water) vs T"
plt.tight_layout()
plt.subplots_adjust( top=0.87 )
plt.savefig( plot_title, dpi=120, transparent=True )

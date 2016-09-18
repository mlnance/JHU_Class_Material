#!/usr/bin/python
__author__="morganlnance"


import matplotlib
matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt
from math import log
import numpy as np


def langmuir( x, Ka ):
    """
    Use the Langmuir adsorption isotherm to calculate X
    :param x: float( adsorption of gas? )
    :param Ka: float( ligand binding constant Kassociation )
    :return float
    """
    X = ( float( Ka ) * float( x ) ) / ( 1 + ( float( Ka ) * float( x ) ) )

    return X


x_val = [ float( 10**-x ) for x in np.arange( 0, 14, 0.1 ) ]
lnx_val = [ log( x ) for x in x_val ]


# vs [x]
plt.plot( x_val, [ langmuir( x, 10**3 ) for x in x_val ], label="Ka = 10^3" )
plt.plot( x_val, [ langmuir( x, 10**6 ) for x in x_val ], label="Ka = 10^6" )
plt.plot( x_val, [ langmuir( x, 10**9 ) for x in x_val ], label="Ka = 10^9" )
plt.plot( x_val, [ langmuir( x, 10**12 ) for x in x_val ], label="Ka = 10^12" )
plt.legend(loc=1)
plt.xlabel( "[x]" )
plt.ylabel( "Langmuir Value" )

# save plot
plot_title = "Langmuir X vs [x] with varying Ka"
plt.tight_layout()
plt.title( plot_title )
plt.subplots_adjust( top=0.87 )
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()


# vs ln[x]
plt.plot( lnx_val, [ langmuir( x, 10**3 ) for x in x_val ], label="Ka = 10^3" )
plt.plot( lnx_val, [ langmuir( x, 10**6 ) for x in x_val ], label="Ka = 10^6" )
plt.plot( lnx_val, [ langmuir( x, 10**9 ) for x in x_val ], label="Ka = 10^9" )
plt.plot( lnx_val, [ langmuir( x, 10**12 ) for x in x_val ], label="Ka = 10^12" )
leg = plt.legend(loc=2, bbox_to_anchor=(1.05,1))
plt.xlabel( "ln[x]" )
plt.ylabel( "Langmuir Value" )

# save plot
plot_title = "Langmuir X vs ln[x] with varying Ka"
plt.tight_layout()
plt.title( plot_title )
plt.subplots_adjust( top=0.87 )
plt.savefig( plot_title, bbox_extra_artists=(leg,), bbox_inches="tight", dpi=120, transparent=True )
plt.close()

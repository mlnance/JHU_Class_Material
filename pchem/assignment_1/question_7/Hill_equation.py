#!/usr/bin/python
__author__="morganlnance"


import matplotlib
matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt
from math import log
import numpy as np


def hill( m, n, pH ):
    """
    Use the modified Hill equation ( m = 1 ) to calculate X
    :param m: float()
    :param n: float()
    :param pH: float( pH )
    :return float
    """
    # the 6 is actually pKa which we are setting to 6 for this question
    X = float(m) * ( (10**(float(n)*(pH - 6))) / (1 + (10**(float(n)*(pH - 6 )))) )

    return X


# pH 0 through 10
pH_vals = np.arange( 0, 10.1, 0.1 )

# m = 1, n = 1, n = 0.75, n = 1.25
plt.plot( pH_vals, [ hill( m=1, n=1, pH=x ) for x in pH_vals ], 'b', label="m=1 n=1" )
plt.plot( pH_vals, [ hill( m=1, n=0.75, pH=x ) for x in pH_vals ], 'r', label="m=1 n=0.75" )
plt.plot( pH_vals, [ hill( m=1, n=1.25, pH=x ) for x in pH_vals ], 'g', label="m=1 n=1.25" )

# n = 1, m = 1, m = 2, m = 3
plt.plot( pH_vals, [ hill( m=1, n=1, pH=x ) for x in pH_vals ], 'b', label="m=1 n=1" )
plt.plot( pH_vals, [ hill( m=2, n=1, pH=x ) for x in pH_vals ], '--r', label="m=2 n=1" )
plt.plot( pH_vals, [ hill( m=3, n=1, pH=x ) for x in pH_vals ], '--g', label="m=3 n=1" )

plt.legend(loc=2)
plt.xlabel( "pH" )
plt.ylabel( "Hill value" )

# save plot
plot_title = "Hill X vs pH with m and n"
plt.tight_layout()
plt.suptitle( plot_title )
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#!/usr/bin/python
__author__="morganlnance"


#import matplotlib
#matplotlib.use( "TKAgg" )
import matplotlib.pyplot as plt
from math import log, factorial, sqrt, pi


def stirlings_approx( N ):
    return ( N * log( N ) ) - N


ii = 2
good_approx = False
diffs = []
x_val = []
while not good_approx:
    x_val.append( ii )
    diff = log( factorial( ii ) ) - stirlings_approx( ii )
    diffs.append( diff )
    error = diff / log( factorial( ii ) ) * 100
    if error < 0.05:
        good_approx = True
        N = ii
    ii += 1
print "N must be about %s for Stirling's approximation to be valid (less than 5%s error).\n" %( N, "%" )

lnsqrt_2piN = [ log(sqrt(2*pi*n)) for n in x_val ]
plt.plot( x_val, diffs, label="Stirling's diff" )
plt.plot( x_val, lnsqrt_2piN, label="ln(sqrt(2*pi*N))" )
plt.xlabel( "N" )
plt.ylabel( "ln( sqrt( 2*pi*N ) )" )
plt.legend( loc=4 )
plot_title = "Difference of Stirling's approximation vs ln( sqrt( 2*pi*N ) )"
plt.tight_layout()
plt.title( plot_title )
plt.subplots_adjust( top=0.87 )
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

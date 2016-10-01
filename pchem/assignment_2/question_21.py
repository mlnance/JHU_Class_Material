#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )

'''
Demonstrates that this is an ideal gas because it roughly obeys Boyles' Law ( PV = m (constant), and P is inversely related to V
At high P, this gas loses ideally, which is to be expected
'''


P = [ 99.49, 101.94, 105.97, 108.44, 116.12, 123.91, 139.06, 149.73, 172.11, 229.4, 91.68, 86.99, 82.7, 77.19, 73.7, 69.4, 66.8, 63.81, 61.39, 58.66, 56.66, 54.1 ]
V = [ 30, 29, 28, 27, 25, 23, 20, 18, 15, 10, 33, 35, 37, 40, 42, 45, 47, 50, 52, 55, 57, 60 ]
inv_V = [ 1.0 / ( float( v ) ) for v in V ]


#plt.plot( V, P )
plt.scatter( V, P )
plt.xlabel( r"V (in$^{3}$)" )
plt.ylabel( "P (Hg)" )
plt.title( "Q21 Pressure vs Volume", fontsize = 16 )
plt.savefig( "Q21_Pressure_vs_Volume", dpi=120, transparent=True )
#plt.show()
plt.close()

#plt.plot( inv_V, P )
plt.scatter( inv_V, P )
plt.xlabel( r"1/V (in$^{3}$)" )
plt.ylabel( "P (Hg)" )
plt.title( "Q21 Pressure vs 1/Volume", fontsize = 16 )
plt.savefig( "Q21_Pressure_vs_Inverse_Volume", dpi=120, transparent=True )
#plt.show()
plt.close()

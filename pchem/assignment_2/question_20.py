#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )

'''
Demonstrates that this is an ideal gas because it obeys Boyles' Law ( PV = m (constant), and P is inversely related to V
'''

P = [ 12, 16, 20, 24, 32, 40, 48 ]  # Hg
V = [ 117.5, 87.2, 70.7, 58.5, 44.2, 35.3, 29.1 ]  # in^3
inv_V = [ 1/v for v in V ]


plt.plot( V, P )
plt.scatter( V, P )
plt.xlabel( r"V (in$^{3}$)" )
plt.ylabel( "P (Hg)" )
plt.title( "Q20 Pressure vs Volume", fontsize = 16 )
plt.savefig( "Q20_Pressure_vs_Volume", dpi=120, transparent=True )
#plt.show()
plt.close()

plt.plot( inv_V, P )
plt.scatter( inv_V, P )
plt.xlabel( r"1/V (in$^{3}$)" )
plt.ylabel( "P (Hg)" )
plt.title( "Q20 Pressure vs 1/Volume", fontsize = 16 )
plt.savefig( "Q20_Pressure_vs_Inverse_Volume", dpi=120, transparent=True )
#plt.show()
plt.close()

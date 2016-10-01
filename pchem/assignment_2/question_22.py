#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 14 } )

'''
This system seems to follow Charles' Law, which is that Volume is proportional to Temperature
This implies that there this is a absolute minimum temperature where motion stops and molcules take up "no" volume
'''


T = [ -124, -100, -75, 0, 25, 30, 45, 100, 225, 323 ]  # Celsius
V = [ 50, 58.05, 66.44, 91.61, 100, 101.68, 106.71, 125.17, 167.11, 200 ]  # mL

plt.plot( T, V )
plt.scatter( T, V )
plt.xlabel( "T (C)" )
plt.ylabel( "V (mL)" )
#plt.xlim( [ 0, 220 ] )
#plt.xticks( range( 0, 221, 20 ) )
plt.title( "Q22 Volume vs Temperature", fontsize = 16 )
plt.savefig( "Q22_Volume_vs_Temperature", dpi=120, transparent=True )
plt.close()


#!/usr/bin/python
__author__="morganlnance"
__question__="ps8_q5"

'''
python random_walk.py
'''


# imports
from random import choice


# step options: 1 step in either direction or nowhere
steps = [ -1, 0, 1 ]

# time to get to a location 100 steps away
location = 0
nsteps = 0
while abs( location ) < 100:
    location += choice( steps )
    nsteps += 1
print "\nIt took %s milliseconds to randomly walk 100 steps in one direction" %( nsteps * ( 10**-3 ) )


# time to get to a location 1000 steps away
location = 0
nsteps = 0
while abs( location ) < 1000:
    location += choice( steps )
    nsteps += 1
print "It took %s milliseconds to randomly walk 1000 steps in one direction\n" %( nsteps * ( 10**-3 ) )

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
tensteps_velocities = []
hundredsteps_velocities = []
while abs( location ) < 100:
    # move a step and count it
    location += choice( steps )
    nsteps += 1
    # calculating velcoity every 10 steps
    if nsteps % 10 == 0:
        tensteps_velocities.append( float( abs( location ) ) / nsteps )
    # calculating velcoity every 100 steps
    if nsteps % 100 == 0:
        hundredsteps_velocities.append( float( abs( location ) ) / nsteps )
print "\n### 100 steps away"
print "It took %s milliseconds to randomly walk 100 steps in one direction" %nsteps
print "   The average velocity using 10-step measurements was %s steps/msec" %round( ( sum( tensteps_velocities ) / len( tensteps_velocities ) ), 5 )
print "   The average velocity using 100-step measurements was %s steps/msec\n" %round( ( sum( hundredsteps_velocities ) / len( hundredsteps_velocities ) ), 5 )



# time to get to a location 1000 steps away
location = 0
nsteps = 0
tensteps_velocities = []
hundredsteps_velocities = []
while abs( location ) < 1000:
    # move a step and count it
    location += choice( steps )
    nsteps += 1
    # calculating velcoity every 10 steps
    if nsteps % 10 == 0:
        tensteps_velocities.append( float( abs( location ) ) / nsteps )
    # calculating velcoity every 100 steps
    if nsteps % 100 == 0:
        hundredsteps_velocities.append( float( abs( location ) ) / nsteps )
print "### 1000 steps away"
print "It took %s milliseconds to randomly walk 1000 steps in one direction" %nsteps
print "   The average velocity using 10-step measurements was %s steps/msec" %round( ( sum( tensteps_velocities ) / len( tensteps_velocities ) ), 5 )
print "   The average velocity using 100-step measurements was %s steps/msec\n" %round( ( sum( hundredsteps_velocities ) / len( hundredsteps_velocities ) ), 5 )

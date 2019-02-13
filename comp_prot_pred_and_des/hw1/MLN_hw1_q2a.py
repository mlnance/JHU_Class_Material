#!/usr/bin/python
__author__="morganlnance"

'''
Homework 1 Question 2a

Usage: ./script.py <nth number of Fibonacci sequence>
Example: ./MLN_hw1_q2a.py 5
'''

###########
# IMPORTS #
###########
import sys


#############
# FUNCTIONS #
#############
def fib(num):
    # num == 0
    # if asking for the 0th fib number, return 0
    if num == 0:
        return 0
    # num == 1
    # if asking for the 1st fib number, return 1
    if num == 1:
        return 1
    # num > 1
    # otherwise, return the sum of fib(num-1) and fib(num-2)
    if num > 1:
        return fib(num-1) + fib(num-2)
    # otherwise, negative number given, return None
    else:
        return None


########
# MAIN #
########
try:
    print "\nFibonacci number %s is %s\n" %(sys.argv[1], fib(int(sys.argv[1])))
except:
    print "\nPlease give me a number for the Fibonacci sequence\n"
    sys.exit()

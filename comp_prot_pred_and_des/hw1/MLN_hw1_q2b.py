#!/usr/bin/python
__author__="morganlnance"

'''
Homework 1 Question 2b

Usage: ./script.py amino.txt
'''

###########
# IMPORTS #
###########
import sys


########
# MAIN #
########
# pull in the amino.txt file from the command line
try:
    amino_file = sys.argv[1]
except IndexError:
    print "\nI need the amino.txt file\n"
    sys.exit()

# read in the amino.txt file
try:
    with open(amino_file, "r") as fh:
        amino_acids = fh.readlines()
except:
    print "\nI couldn't open and/or read your argument\n"
    sys.exit()

# print the first three letters of each amino acid
for aa in amino_acids:
    print aa[:3].strip()

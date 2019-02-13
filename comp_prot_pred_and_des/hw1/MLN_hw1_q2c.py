#!/usr/bin/python
__author__="morganlnance"

'''
Homework 1 Question 2c

Usage: ./script.py
'''

###########
# IMPORTS #
###########
import sys


###########
# CLASSES #
###########
class Person:
    def __init__(self, name, phone_num, email):
        self.name = name
        self.phone_num = phone_num
        self.email = email


class Contacts:
    # can initialize with or without a starting Person contact
    def __init__(self, person_ = None):
        self.contacts = []
        if person_ is not None:
            self.contacts.append(person_)

    def add_contact(self, person_):
        self.contacts.append(person_)

    def add_contacts(self, persons_):
        for person_ in persons_:
            self.contacts.append(person_)

    def print_contacts(self):
        for p in self.contacts:
            print "%s\t%s\t%s" %(p.name, p.phone_num, p.email)

    def sort_contacts(self):
        # assumes p.name is of type "First Last"
        last_names = [p.name.split()[-1] for p in self.contacts]
        last_names.sort()
        local_contacts = []
        for last_name in last_names:
            for p in self.contacts:
                if last_name == p.name.split()[-1]:
                    print "%s\t%s\t%s" %(p.name, p.phone_num, p.email)
                    local_contacts.append(p)
                    break
        self.contacts = [ii for ii in local_contacts]

    def phone_lookup(self, name):
        for person in self.contacts:
            if person.name.lower() == name.lower():
                print "%s --> %s" %(person.name, person.phone_num)
                return True
        print "\n%s does not exist in Contacts\n" %name
        return None

    def email_lookup(self, name):
        for person in self.contacts:
            if person.name.lower() == name.lower():
                print "%s --> %s" %(person.name, person.email)
                return True
        print "\n%s does not exist in Contacts\n" %name
        return None


########
# MAIN #
########
contact_manager = Contacts()
p1 = Person("Jane Richardson", "123-555-7891", "jprotein@duke.edu")
p2 = Person("Cyrus Levintal", "555-555-8910", "paradox@gmail.com")
p3 = Person("Shoshana Wodak", "819-555-1234", "docking@sickkids.edu")
p4 = Person("David Baker", "425-555-3341", "hiker@uw.edu")
contact_manager.add_contacts([p1, p2, p3, p4])

print "\n1) Example of sorting"
print "Starting contacts list:"
contact_manager.print_contacts()
print "\nAfter sorting:"
contact_manager.sort_contacts()
print "\n\n2) Example of phone lookup by name"
contact_manager.phone_lookup("jane richardson")
print "\n\n3) Example of email lookup by name"
contact_manager.email_lookup("david baker")
print

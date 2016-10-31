#####################################################################
#
# Colour by object
#
#####################################################################
 
def set_torsion1():
 
        """
        
        this function sets bond of two named atoms to be moved by torsion command

        """
        cmd.edit("x1","y1")

def set_torsion2():

        """

        this function sets bond of two named atoms to be moved by torsion command

        """
        cmd.edit("x2","y2")

def set_torsion3():

        """

        this function sets bond of two named atoms to be moved by torsion command

        """
        cmd.edit("x3","y3")

def set_torsion4():

        """

        this function sets bond of two named atoms to be moved by torsion command

        """
        cmd.edit("x4","y4")

def set_torsion5():

        """

        this function sets bond of two named atoms to be moved by torsion command

        """
        cmd.edit("x5","y5")


# Add to the PyMOL command list 
cmd.extend("set_torsion1",set_torsion1)
cmd.extend("set_torsion2",set_torsion2)
cmd.extend("set_torsion3",set_torsion3)
cmd.extend("set_torsion4",set_torsion4)
cmd.extend("set_torsion5",set_torsion5)


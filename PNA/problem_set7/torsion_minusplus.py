#####################################################################
#
# make new commands to rotate bonds in plus and minus direction
#
#####################################################################
 
def torsion_plus():
 
        """
        
        this function moves the selected bond (with edit and pk1 pk2) 
        in a positive direction. 
        This is a hack because I keep getting the error:

        Traceback (most recent call last):
        File "/users/warren/pymol/products/MacPyMOL.app/pymol/modules/pymol/parser.py", 
        line 188, in parse
        result=apply(kw[nest][0],args[nest],kw_args[nest])
        File "/users/warren/pymol/products/MacPyMOL.app/pymol/modules/pymol/cmd.py", 
        line 808, in _special
        apply(my_special[k][1],my_special[k][2],my_special[k][3])
        TypeError: torsion() takes exactly 1 argument (2 given)bond to be moved by torsion command

        """
        cmd.torsion("10")

def torsion_minus():

        """

        this function moves the selected bond (with edit and pk1 pk2)         
        in a positive direction.
        This is a hack because I keep getting the error:

        Traceback (most recent call last):
        File "/users/warren/pymol/products/MacPyMOL.app/pymol/modules/pymol/parser.py",         line 188, in parse
        result=apply(kw[nest][0],args[nest],kw_args[nest])
        File "/users/warren/pymol/products/MacPyMOL.app/pymol/modules/pymol/cmd.py",         line 808, in _special
        apply(my_special[k][1],my_special[k][2],my_special[k][3])
        TypeError: torsion() takes exactly 1 argument (2 given)bond to be moved by torsion command
        """

        cmd.torsion("-10")


# Add to the PyMOL command list 
cmd.extend("torsion_plus",torsion_plus)
cmd.extend("torsion_minus",torsion_minus)


#!/usr/bin/python
__website__="http://azevedolab.net/resources/dihedral_angle.pdf"


def initial_vectors(p1_in, p2_in, p3_in, p4_in):
    """Function to set up initial vectors"""
    import numpy as np
    # Set initial values for arrays
    p1 = np.zeros(3)
    p2 = np.zeros(3)
    p3 = np.zeros(3)
    p4 = np.zeros(3)
    # Set initial coordinates (http://www.stem2.org/je/proteina.pdf)
    p1[:] = p1_in
    p2[:] = p2_in
    p3[:] = p3_in
    p4[:] = p4_in
    return p1,p2,p3,p4 


def calc_q_vectors(p1,p2,p3,p4):
    """Function to calculate q vectors"""
    import numpy as np
    # Calculate coordinates for vectors q1, q2 and q3
    q1 = np.subtract(p2,p1) # b - a
    q2 = np.subtract(p3,p2) # c - b
    q3 = np.subtract(p4,p3) # d - c
    return q1,q2,q3


def calc_cross_vectors(q1,q2,q3):
    """Function to calculate cross vectors"""
    import numpy as np
    # Calculate cross vectors
    q1_x_q2 = np.cross(q1,q2)
    q2_x_q3 = np.cross(q2,q3)
    return q1_x_q2, q2_x_q3


def calc_nornals(q1_x_q2,q2_x_q3):
    """Function to calculate normal vectors to planes"""
    import numpy as np
    # Calculate normal vectors
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))
    return n1,n2


def calc_orthogonal_unit_vectors(n2,q2):
    """Function to calculate orthogonal unit vectors"""
    import numpy as np
    # Calculate unit vectors
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2,q2)))
    u2 = np.cross(u3,u1)
    return u1,u2,u3


def dihedral_angle(n1,u1,u2,u3):
    """Function to calculate dihedral angle"""
    import numpy as np
    import math
    # Calculate cosine and sine
    cos_theta = np.dot(n1,u1)
    sin_theta = np.dot(n1,u2)
    # Calculate theta
    theta = -math.atan2(sin_theta,cos_theta) # it is different from Fortran math.atan2(y,x)
    theta_deg = np.degrees(theta)
    # Show results
    print("theta (rad) = %8.3f"%theta)
    print("theta (deg) = %8.3f"%theta_deg)



def calc_dihedral_angle(p1_in, p2_in, p3_in, p4_in):
    """
    Calculate the dihedral angle between four points in space
    p*_in = [x, y, z]
    """
    # Call initial_vectors() functions
    p1,p2,p3,p4 = initial_vectors(p1_in, p2_in, p3_in, p4_in)
    # Call calc_q_vectors(p1,p2,p3,p4) function
    q1,q2,q3 = calc_q_vectors(p1,p2,p3,p4)
    # Call calc_cross_vectors(q1,q2,q3) function
    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1,q2,q3)
    # Call calc_nornalss(q1_x_q2,q2_x_q3) function
    n1, n2 = calc_nornals(q1_x_q2,q2_x_q3)
    # Call calc_orthogonal_unit_vectors(n2,q2) function
    u1,u2,u3 = calc_orthogonal_unit_vectors(n2,q2)
    # Call calc_dihedral_angle(u1,u2,u3) function
    dihedral_angle(n1,u1,u2,u3)


#if __name__ == "__main__":
#    calc_dihedral_angle()

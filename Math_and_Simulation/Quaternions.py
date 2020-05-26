#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:05:55 2018

@author: top40ub
"""

"""
In Navigation in 3D is a highly problematic and non trivial task. Due to the 
anisotropic definition von spherical coordinate systems rotations and movement
along isolines or fixed points is tricky. A way to do this a bit easier is to
an old and archaric concept from paricle physics quaternions. To be true their
mathematical formulation excites those of vectors in therms of age.
"""


import numpy as np


"""
Transformation of 3d vectors to 4 dimensional quaternions
Input : np.array([x,y,z])
Outout : np,array([0,x,y,z])
"""
"""
Function name : vec_to_quat()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def vec_to_quat(vec):
    if vec.shape[0] == 3:
        quat = np.append([0], vec)
        return quat

    else:
        raise Exception('Error the vector is not a 3D vector [x,y,z]')

    	
"""
Calculate the conjugate quaternion of given quaternion
Input : np.array([w,x,y,z])
Outout : np,array([w,-x,-y,-z])
"""
"""
Function name : quat_conjug()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def quat_conjug(quat):
    if quat.shape[0] == 4:
        w, x, y, z = quat
        quatconjung = np.array([w, -x, -y, -z])
        return quatconjung
    else:
        raise Exception('Error the quaternion is not a 4D [w,x,y,z]')


"""
Normalisation of given vector with a certain tolerance
Input : np.array([x,y,z])
Outout : np,array([x,y,z]) with x² + y² + z² = 1
"""
"""
Function name : vec_normal()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def vec_normal(vec, tolerance=0.00001):
    if vec.shape[0] == 3:
    	vec_len_quad = np.sum([n * n for n in vec])
    	if abs(vec_len_quad - 1.0) > tolerance:
        	vec_len = np.sqrt(vec_len_quad)
        	vec = vec/vec_len
    	return vec
    else:
        raise Exception('Error the vector is not a 3D vector [x,y,z]')


"""
This function defins how to quaternions multiplicate. The multiplication core
can be interpreted as a kind of rotation. For deeper information read the
coresponding wikipedia page or one of the plenty graphical visualition books
about computer vision like I' ve done :)
Input : np.array([w1,x1,y1,z1]), np,array([w2,x2,y2,z2])
Outout : np,array([w_mult,x_mult,y_mult,z_mult])
"""
"""
Function name : quat_mult()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def quat_mult(quat1, quat2):
    if quat1.shape[0] == 4 and quat2.shape[0] ==4:
        w1, x1, y1, z1 = quat1
        w2, x2, y2, z2 = quat2
        w_mult = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x_mult = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y_mult = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
        z_mult = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
        quat_mult = np.array([w_mult, x_mult, y_mult, z_mult])
        return quat_mult
    else:
        raise Exception('Error one quaternion is not a 4D [w,x,y,z]')
        

"""
Quaternionic multiplication between quaternion and vector. The vector
is transformed. A second quaternionc multiplication with the conjugate of given
quaternion returns the three vector compoments [1:] of the operation.
Input : np.array([x1,y1,z1]), np,array([w2,x2,y2,z2])
Outout : np,array([w_mult,x_mult,y_mult,z_mult])
"""
"""
Function name : qaut_vec_mult()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def qaut_vec_mult(vec1, quat1):
    if vec1.shape[0] == 3 and quat1.shape[0] ==4:
    	quat2 = vec_to_quat(vec1)
    	return quat_mult(quat_mult(quat1, quat2), quat_conjug(quat1))[1:]
    else:
        raise Exception('Error the quaternion is not a 4D [w,x,y,z] or the vector is not 3D [x,y,z]')
        

"""
Calcuation to transform a unit vector or any vector that will be normalsied 
into a quaternion that elements are object to a rotation around half theta.
!Theta needs to be of type degree!
Input : np.array([x,y,z]), np.array([theta])
Outout : np,array([w_a,x_a,y_a,z_a])
"""
"""
Function name : axisangle_to_quat()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def axisangle_to_quat(vec, theta):
    if vec.shape[0] == 3 and theta.shape[0] == 1:
        vec = vec_normal(vec)
        x, y, z = vec
        theta = theta[0]/2
        theta = np.deg2rad(theta)
        w_a = np.cos(theta)
        x_a = x * np.sin(theta)
        y_a = y * np.sin(theta)
        z_a = z * np.sin(theta)
        quat_a = np.array([w_a, x_a, y_a, z_a])
        return quat_a
    else:
        raise Exception('Error the vector is not 3D [x,y,z] or the angle 1D [theta]')
        
        
"""
The opposite function of axisangle_to_quat, it returns the axisangle the 
corresponding unit vector of a given quaternion.
Input : np.array([w,x,y,z])
Outout : np,array([x_a,y_a,z_a]), np.arra([theta])
"""
"""
Function name : quat_to_axisangle()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def quat_to_axisangle(quat):
   if quat.shape[0] == 4: 
    	w, vec = quat[0], quat[1:]
    	theta = np.arccos(w) * 2.0
    	return vec_normal(vec), theta
    
   else:
       raise Exception('Error the quaternion is not a 4D [w,x,y,z]')


"""
Main function for the quaternionic rotation of 3D vectors:
A vector (vec) is rotated around a second vector a_unit with angle theta. 
a_unit is a axisvector but mostly a basis vector of given cartesic coordinate
system. In case of a cartesic basis vector we have the rotation around 
coordinate axsis.
Input : np.array([x,y,z]), np.array([theta]), 
        np.array(x2,y2,z2) with x2² + y2² + z2² = 1
Output : np.array([x1_rot,y1_rot,z1_rot])
"""
"""
Function name :vec_angle_rot()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a)...
	b)...
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def vec_angle_rot(vec,theta,a_unit):
    if vec.shape[0] == 3 and theta.shape[0] == 1 and vec.shape[0] == 3:
    	quat1 = vec_to_quat(vec)
    	quat_a = axisangle_to_quat(a_unit,theta)
    	v_rot = quat_mult(quat_mult(quat_a,quat1),quat_conjug(quat_a))[1:]
    	return v_rot
    else:
        raise Exception('Error one vector is not 3D [x,y,z] or the angle 1D [theta]')


if __name__=='__main__':

    vector = np.array([[0,0,1],
                       [1,0,1],
                       [-1,0,1],
                       [0,1,1],
                       [0,-1,1]])
    print(np.linalg.norm(vector,axis = 1))
    vector = (vector.T/np.linalg.norm(vector,axis = 1)).T
    xu = np.array([0,1,0])
    zu = np.array([0,0,1])

    
    pol = np.array([90])
    az = np.array([45])
    for i,j in zip(vector,range(vector.shape[0])):
        i = vec_angle_rot(i,pol,xu)
        i = vec_angle_rot(i,az,zu)
        vector[j] = i
        print(i)
    #vector = np.round(vector)
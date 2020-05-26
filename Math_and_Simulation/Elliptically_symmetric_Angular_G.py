#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 21:30:59 2019

@author: top40ub
"""

import numpy as np
from scipy.stats import norm
from numpy import pi, cos, sin, arccos, arange
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from time import time as time_it

import Math_and_Simulation.Quaternions as Quaternions
from Math_and_Simulation import Pdf_Convolver as PdfC


"""
Function name : sphere_vector()
***Description***

Inspired by the field of spherical distributions and navigation I create a
distribution an the 3d unit sphere. The distribution kernel for this special 
distribution is calcualted from a 3D Multivariate Normal Gaussian Distribution.
Imput parameters are: 
    aperture: see build_sphereical_grid
    angle : see build_spherical_grid
    CV : Covariance, the Metric or the kernel for the distribution
    mµ : mean vector of the distribution
    montecarlo_iter : maximum number of tries to find a new direction
    offset = see build_spherical_grid
    grid_point = see build_spherical_grid
    plot = Boolean decision for a 3D plot of the spherical distribution
The function returns an unit vector, np.array = [x,y,z]

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
def sphere_vector(aperture, angle, CV, mµ, montecarlo_iter, offset = 0, points_on_sphere = 4000, plot = False):
    "new seed for random decisions and RNG"
    
    np.random.seed()
    grid_points, sampling_factor = sample_gridpoint_number(aperture, offset, points_on_sphere)
    vector_l = build_spherical_grid(aperture, angle, offset, grid_points)
    
    fag_norm = build_fAG_norm(vector_l,CV,mµ)
    alpha = build_alpha(vector_l,CV,mµ)
    exp_part = build_exp_part(alpha,CV,mµ)
    m2_norm = build_M2_norm(alpha)

    """
    This is the spherical distribution (fAG) with its probability 
    calculated at each grid node of the grid map
    The function returns an array, np.array = [fAG0,..,fAGn]
    """
    fAG = fag_norm * exp_part * m2_norm / sampling_factor
    sum_fAG= np.sum(fAG)
    fAG = fAG / sum_fAG
    """
    This part of the function runs a Monte Carlo simulation to draw a random
    vector out of the sperical distribution.
    """
    
    a = monte_carlo(fAG, grid_points, montecarlo_iter)
    if plot == True:
        plot_sphere_vector(vector_l,a,fAG)

    return vector_l[a]

"""
Function name : sample_gridpoint_number()
***Description***

--- plain text ---

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
def sample_gridpoint_number(aperture,offset, points_on_sphere):
    r = 1
    A = 2 * np.pi * r * r * (1 - np.cos(np.deg2rad(aperture+offset)))
    A = A - 2 * np.pi * r * r * (1 - np.cos(np.deg2rad(offset)))
    grid_sampling =points_on_sphere/(4* np.pi * r * r)
    grid_points =  int(A * grid_sampling)
    return grid_points, grid_sampling


"""
Function name : build_spherical_grid()
***Description***

This function creates a nearly equidistant grid on the 3d unit sphere. The used
method distributes grid points starting from the north pole in Fibonacci series
to obtain a full coverage. The grid resembles entangled Fibonacci series like
one see in sunflowers. The paramters are:
    gridpoints : number of grid nodes [0,..,n]
    offset : instead of the northpol the series start at latitude = offset
    aperture : the nodes are within latitude [offset, offset + aperture]
    angle : quaternionic rotate all grid point in a way that the new north pole
            is placed at [polar, azimutal] coordinate angle = [pol,az]
With the quaternionic rotation the suface area of a solid angle can be placed 
at each position.
The function returns the nodes as np.array = [[x0,y0,z0],...,[xn,yn,zn]]

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
def build_spherical_grid(aperture, angle, offset, grid_points):    
    pol = np.array([angle[0]])
    az = np.array([angle[1]])
    indices = arange(0, grid_points, dtype=float) + 0.5
    aperture_vol = (1-cos(np.deg2rad(aperture)))/2
    phi =np.deg2rad(offset) + arccos(1 - 2*indices/grid_points*aperture_vol)
    theta = pi  * indices * (1 + 5**0.5)
    vector_l = np.array([cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)])
    xu = np.array([0,1,0])
    zu = np.array([0,0,1])
    vector_l = vector_l.T
    
    for i in range(grid_points):
        vector_l[i] = Quaternions.vec_angle_rot(vector_l[i],pol,xu)
        vector_l[i] = Quaternions.vec_angle_rot(vector_l[i],az,zu)
    unit_vector_array = vector_l
    return unit_vector_array


"""
Function name : build_fAG_norm()
***Description***

This function constructs the normalization part of the spherical distribution.
An array with the grid node coordinates (vector_l) a Covariance Matirx (CV) and 
a mean vector (mµ) are the parameters. The returned objct is a array for the 
normalisation at each node. 
The function returns an array, np.array = [fag_norm0,..,fag_norm1]

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
def build_fAG_norm(vector_l, CV, mµ):
    cd_norm = 1/(2*np.pi)/np.sqrt(np.linalg.det(CV))
    inv_CV = np.linalg.inv(CV)
    fag_norm = np.array([cd_norm * 1/np.dot(np.dot(t.T,inv_CV),t)**(3/2)    for t in vector_l])
    return fag_norm



"""
Function name : build_alpha()
***Description***

This function calculates the value alpha for each grid node in the grid
node postion array (vector_l) out of the covarince (CV) and 
the mean vector (mµ).
The function returns an array, np.array = [alpha0,..,alphan] 
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
def build_alpha(vector_l,CV,mµ):
    inv_CV = np.linalg.inv(CV)
    alpha = np.array([np.dot(np.dot(t.T,inv_CV),mµ)/np.dot(np.dot(t.T,inv_CV),t)**(1/2)   for t in vector_l])
    return alpha


"""
Function name : build_exp_part()
***Description***

With alpha, the covarinace (CV) and the mean vector (mµ) the exponential part
of the distribution is calculated. Each entry in the return array corresponds
to one node in the node_position array (vector_l).
The function returns an array, np.array = [exp_part0,..,exp_partn]

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
def build_exp_part(alpha,CV,mµ):
    inv_CV = np.linalg.inv(CV)
    exp_part = np.exp(alpha**2 - np.dot(np.dot(mµ.T,inv_CV),mµ))
    return exp_part


"""
Function name : build_M2_norm()
***Description***

Since the spherical distribution is a marginal distribution of the full 3D
Multivarite Normal Gaussian Distribution the integration over the the radial
component is not trivial in case of non zero mean (mµ). Therfore a second
normalisation is performed that takes the distribution of the alpha values into
account
The function returns an array, np.array = [m2_norm0,..,m2_normn]

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
def build_M2_norm(alpha):
    m2_norm = (1+alpha**2)*norm.cdf(alpha) + alpha*norm.pdf(alpha)
    return m2_norm


"""
Function name : plot_sphere_vector()
***Description***

This function plot as a spherical distribution (fAG) at the grid nodes in grid
node position arrray (vector_l). a is the the position of the node that was
drawn over he distribution via a Monte Carlo Simultation.
The plot draws a the 3D plot with: 
    Grid point colour indicating their probability.
    The Cartesic coordinate axes (x,y,z) in green, blue and red
    The Monte Carlo Simulated unit vector at a in magenta
The function creates a 3D plot

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
def plot_sphere_vector(vector_l, a, fAG):
        fig = plt.figure()
        [x, y, z] = vector_l.T

        listecolor =np.array([c for i,j,k,c in zip(x,y,z,fAG)])
        #listecolor = fAG.reshape(grid_points)
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.quiver(1, 0, 0,.5,0,0, color = 'g')
        ax1.quiver(0, 1, 0,0,.5,0, color = 'b')
        ax1.quiver(0, 0, 1,0,0,.5, color = 'r')
        ax1.quiver(*vector_l[a],*vector_l[a]*2, color = 'm')
        ax1.scatter(x,y,z,c = listecolor)
        ax1.set_xlim(-1,1)
        ax1.set_ylim(-1,1)
        ax1.set_zlim(-1,1)

        ax1.set_xlabel('$X$', fontsize=10)
        ax1.set_ylabel('$Y$', fontsize=10)
        ax1.set_zlabel('$Z$', fontsize=10)
        ax1.set_xticks([-1,0,1])
        ax1.set_yticks([-1,0,1])
        ax1.set_zticks([-1,0,1])
        ax1.view_init(elev=-3,azim=-75+180)
        plt.show()



"""
Function name : rotation_mat()
***Description***

This is a function the construct a rotation matrix for the three rotation needed
in the Euler rotation picture.
These matrices are used for the covarinace construction
Theta must be aa array with the rotation angles in degree np.array([z,x,z])
The function returns a matrix, np.array([[a,c,b],[d,e,f],[g,h,i]])

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
def rotation_mat(theta):
    """
    To stay in the Euler rotation picture we rotate: 
    First around Rz this is the azimutal angle tz: (0,2pi) or (0,2pi) 
    // counter clockwise looked from -z direktion
    Then around Rx' which is the polar angle  tx: (0,pi) 
    After this rotation the z' axsis in the new coordinate frame is fixed. 
    // counter clockwise looked from -x direktion
    The last rotation is again around the a z axsis, the new z' axsis. 
    This will put the frame x and y axes in the right orientation. 
    // counter clockwise looked from -z direktion
    """
    tz,tx,tz2 = np.deg2rad(theta)
 
    Rz = np.array([[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0,0,1]])
    Rx = np.array([[1,0,0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
    Rz2 = np.array([[np.cos(tz2), -np.sin(tz2), 0], [np.sin(tz2), np.cos(tz2), 0], [0,0,1]])
    
    'return the combiend rotation matrix'
    return np.dot(np.dot(Rz2,Rx), Rz)


"""
Function name : print_pos_vec()
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
def print_pos_vec(vector_l,fAG, pos_vec):
    fig = plt.figure()
    [x, y, z] = vector_l.T #*np.array([[5,3,2]]).T

    listecolor =np.array([c for i,j,k,c in zip(x,y,z,fAG)])
    listecolor = fAG.reshape(grid_points)
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.quiver(1, 0, 0,.5,0,0, color = 'g')
    ax1.quiver(0, 1, 0,0,.5,0, color = 'b')
    ax1.quiver(0, 0, 1,0,0,.5, color = 'r')
    ax1.quiver(*pos_vec,*pos_vec*2, color = 'm')
    ax1.scatter(x,y,z,c = listecolor)
    ax1.set_xlim(-1,1)
    ax1.set_ylim(-1,1)
    ax1.set_zlim(-1,1)

    ax1.set_xlabel('$X$', fontsize=10)
    ax1.set_ylabel('$Y$', fontsize=10)
    ax1.set_zlabel('$Z$', fontsize=10)
    ax1.set_xticks([-1,0,1])
    ax1.set_yticks([-1,0,1])
    ax1.set_zticks([-1,0,1])
    ax1.view_init(elev=-3,azim=-75+180)
    plt.show()


def monte_carlo(fAG, grid_points, montecarlo_iter):
    for i in range(int(montecarlo_iter)):

        a = random.randrange(grid_points)
        rand_prob= np.random.random() 
        calc_prob= fAG[a]
        if (rand_prob <= calc_prob):
            return a
    print('FAIL',  np.amax(fAG))
 
def monte_carlo2(fAG, grid_points):
    a = np.random.choice(np.arange(grid_points),p = fAG)

    return a


"""
Function name : start_point_generator()
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
def start_point_generator(start_pos_list,
                          start_angle_list,
                          aperture = 5, offset = 0,
                          points_on_sphere = 1000, 
                          metric = np.array([[0,0,0]]),
                          common_starting_node = True,
                          tangent_orientation = 'none'):
    start_node_counter = 1
    start_pos = np.array([])
    start_angle = np.array([])
    common_node = np.array([])
    for p,a in zip(start_pos_list,start_angle_list):

        node_number, sampling = sample_gridpoint_number(aperture, offset, points_on_sphere)
        pos = build_spherical_grid(aperture, a ,offset, node_number)

        directions = np.array([PdfC.convert_vector_angle(vec)[1][0] for vec in pos])
        if tangent_orientation== 'pol':

            vec_pol= np.array([np.cos(np.deg2rad(directions[:,0]))*np.cos(np.deg2rad(directions[:,1])),
                             np.cos(np.deg2rad(directions[:,0]))*np.sin(np.deg2rad(directions[:,1])),
                             -np.sin(np.deg2rad(directions[:,0]))]).T
            directions = np.array([PdfC.convert_vector_angle(vec)[1][0] for vec in vec_pol])
    
        if tangent_orientation == 'azi':

            vec_az= np.array([-np.sin(np.deg2rad(directions[:,1])),
                       np.cos(np.deg2rad(directions[:,1])),
                       0*np.deg2rad(directions[:,1])]).T
            directions = np.array([PdfC.convert_vector_angle(vec)[1][0] for vec in vec_az])
            
        node_counter = np.ones(node_number) * start_node_counter
        pos = pos*metric + p

        start_pos = np.append(start_pos, pos)

        start_angle = np.append(start_angle,directions)
 
        if common_starting_node == True:
            start_node_counter += 1
            common_node = np.append(common_node,node_counter)
        else:
            node_counter = np.arange(start_node_counter,pos.shape[0]+start_node_counter)
            common_node = np.append(common_node,node_counter)
            start_node_counter = pos.shape[0]+start_node_counter
    start_pos = start_pos.reshape(int(start_pos.shape[0]/3),3)
    start_angle = start_angle.reshape(int(start_angle.shape[0]/2),2)
    common_node = np.array([common_node]).T
    
    return start_pos, start_angle, common_node

if __name__ == '__main__':
    time = 0
    points_on_sphere = 2000
    aperture = 20
    offset =65
    angle=np.array([0, 0])
    
    
    cv = np.array([[0.1,0,0],[0,5,0],[0,0,5]])
    rot = rotation_mat(np.array([0,0,0]))
    cv = np.dot(np.dot(rot,cv),rot.T)
    print(cv)
    eig = np.linalg.eigh(cv)
    cl = (eig[0][2] - eig[0][1])/np.sum(eig[0])
    cp = 2*(eig[0][1] - eig[0][0])/np.sum(eig[0])
    cs = 3*eig[0][0]/np.sum(eig[0])
    ca = cl + cp

    mµ = np.array([0,0,0])+ np.array([0,0,0])+np.array([0,0,0])
    pos_vec = np.array([0,0,0])
    for i in range(1):
        tic = time_it()

        grid_points, sampling_factor = sample_gridpoint_number(aperture,offset, points_on_sphere)
        vector_l = build_spherical_grid(aperture, angle, offset, grid_points)
        vector_l_s = a = vector_l * np.array([1,1,1])
        vector_l = vector_l_s/np.array([np.linalg.norm(vector_l_s,axis = 1)]).T
        fag_norm = build_fAG_norm(vector_l,cv,mµ)
        alpha = build_alpha(vector_l,cv,mµ)
        exp_part = build_exp_part(alpha,cv,mµ)
        m2_norm = build_M2_norm(alpha)
        fAG = fag_norm * exp_part * m2_norm
        fAG_prob = fAG / sampling_factor
        fAG_prob = fAG_prob/ np.sum(fAG_prob)
    
        
        a = monte_carlo2(fAG_prob, grid_points)
        pos_vec = pos_vec+ vector_l[a]
        toc = time_it()
        time += toc-tic
    print('time:', time, 'in sec')
    pos_vec = pos_vec/np.linalg.norm(pos_vec)
    print(pos_vec)
    print_pos_vec(vector_l,fAG_prob, pos_vec)
    pos_vec = np.array([0,0,0])
    montecarlo_iter = 10000
    time2 = 0
    for i in range(1):
        tic = time_it()
        pos_vec = sphere_vector(aperture, angle, cv, mµ, montecarlo_iter, offset = offset, points_on_sphere = 1000, plot = False)
        toc = time_it()
        time2 += toc-tic
    print('time2:', time2, 'in sec')
    pos_vec = pos_vec/np.linalg.norm(pos_vec)
    print(pos_vec)
    print_pos_vec(vector_l,fAG_prob, pos_vec)
    

    start_pos_list = np.array([[125,125,25],
                               [125,125,225],
                               [125,25,125],
                               [125,225,125],
                               [25,125,125],
                               [225,125,125]])
   
    start_angle_list = np.array([[0,0],
                                 [0,0],
                                 [90,90],
                                 [90,90],
                                 [90,0],
                                 [90,0]])  
    start_pos, start_angle,\
    common_node = start_point_generator(start_pos_list,
                                        start_angle_list,
                                        aperture = 10, offset = 90,
                                        points_on_sphere = 200, 
                                        metric = np.array([[1,1,1]]),
                                        common_starting_node = True,
                                        tangent_orientation = 'azi')

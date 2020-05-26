#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:05:55 2018

@author: top40ub
"""

import numpy as np

'This python file runs a monte carlo process over a given pdf and aperture/solid angle to draw the next direction'
import Math_and_Simulation.Elliptically_symmetric_Angular_G as fAG

import Math_and_Simulation.Quaternions as Quaternions


"""
The function pdf_convolver calculates with the force_drift and force_struc 
together with the internal parameter from cod_object the covariance and the mean
for the distribution kernel of the sphere_vector function
The paramters are:
    cod_obj : a GrowthCone Class object
    force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 componetes
                  for the diffusion covarince

To caculate a new direction the kernel and mean of a pdf 
probability density function is constructed out of sum of 3 parts
The pdf is a multivariate normal distribuition with covariance SIGMA and mean µ
    a) µ = (0,0,0) . symmetric along each eigenvector of its SIGMA 
       SIGMA  = covariance : out of force_struc the structure tensor
    b) µ = drift vector : out of the force_drift  
       SIGMA = Diag(1,1,1)*sigma = 0.01 : sharp discret signal, nearly a delta peak
    c) µ = cod.pols_list[-1] : internal parameter, last growthdirection 
       SIGMA = diag(1,1,1)*0.2 : internal parameter, Correlation of growth diection with last directions
       
The function returns a normalised unit_vector new_posvec, np.array([x,y,z])
"""
"""
Function name : pdf_convolver()
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
def pdf_convolver(cod_obj, force_drift, force_struc):

    'Drift vector components'
    'Drift signal//all directions'
    driftx, drifty, driftz =  force_drift
    
    'Drift mean µ-vector'
    drift_mu = np.array([driftx,drifty,driftz])
    if cod_obj.max_drift != 0: 
        drift_mu = drift_mu/cod_obj.max_drift
    #drift_mu = - drift_mu #to move away from borders instead getting attracted

    
    'Drift/shift covariance //SIGMA'
    struc_ten_drift = np.array([[1,0,0],[0,1,0],[0,0,1]]) * 0.01

    'Structure tensor components'
    'Diffusion metric // Diffusive tensor components'
    diffxx = force_struc[0]
    diffxy = force_struc[1]
    diffxz = force_struc[2]
    diffyy = force_struc[3]
    diffyz = force_struc[4]
    diffzz = force_struc[5]
    
    'Diffusion mena µ-vector'
    diff_mu = np.array([0,0,0])
    
    'Diffuison covariance //SIGMA'    
    struc_ten_diff = np.array([[diffxx,diffxy,diffxz],[diffxy,diffyy,diffyz],[diffxz,diffyz,diffzz]])
    'Inverted length scales for real space opjects'
    eigv_diff = struc_ten_eig(struc_ten_diff,cod_obj.max_eig, cod_obj.min_eig, cod_obj.proxy_reverse_eig)
    struc_ten_diff = revers_eig(eigv_diff)

    'Correlation mean µ-vector'
    corr_mu = cod_obj.vec_mem
    
    'Correlation covariance //SIGMA'
    struc_ten_corr = np.array([[1,0,0],[0,1,0],[0,0,1]])*0.01
    
    """Proxy pdf parameter"""
    if cod_obj.proxy_drift == True:
        drift_mu = np.array([0,0,0]) #Proxy drift mean vector drift_mu
    if cod_obj.proxy_tensor == True:
        struc_ten_diff= np.array([[0.5,0,0],[0,0.5,0],[0,0,1]]) #Proxy diffusion tensor
    if cod_obj.proxy_corr == True:
        corr_mu = np.array([0,0,0]) #Proxy correlation mean vector corr_mu
    """Calculate all sums:"""
    
    'Sum of all three covariances'
    struc_ten_sum = struc_ten_drift +struc_ten_diff + struc_ten_corr

    'Sum of all three means' 
    shift_mu = drift_mu + diff_mu + corr_mu
    """
    A Monte Carlo Simulation draws the new direction from a distribution 
    constructed out of the covarinaces and means. The cones internal parameter
    like the aperture restricts the new growth direction inside a solid angle.
    """
    new_posvec = fAG.sphere_vector(cod_obj.aperture, 
                                   cod_obj.angle_list[-1], 
                                   struc_ten_sum, shift_mu, 
                                   cod_obj.montecarlo_iterations, 
                                   offset = 0, 
                                   points_on_sphere = 2000, 
                                   plot = False)

    return new_posvec/np.linalg.norm(new_posvec)


"""
This function takes the covarinace matirx that is the sum of all three parts
drift, diffusion and correlation and performs a diagonalisation. This results
in three eigenvalues sorted by magnitude.
The growth algorithm uses the covarinace as structual information of the substrate
enhanced by the internal growth parameter of peristence and correlation and by
the external force like the drift. Therfore the eigenvalues musst be inverted
in respect to the other eigenvalues. In 2D is is easily done with transformation
matrix [[0,-1],[1,0]] but in 3D one has to take the quotient of the middle 
eigenvalue towards the other into account.

The returned object is a tuple the same way as function np.linalg.eigh
"""
"""
Function name : struc_ten_eig()
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
def struc_ten_eig(struc_ten, max_eig, min_eig, proxy_reverse_eig):
    eig = np.linalg.eigh(struc_ten)
    epsilon = 0.1

    for i in range(3):
        if eig[0][i] < epsilon:
            eig[0][i]= epsilon
    if proxy_reverse_eig == True:
        sum_eig =1/(1+np.sum(eig[0])) #np.linalg.norm(eig[0])
        ten1 = (eig[0][0]/sum_eig)**(-1)*(eig[0][0]/sum_eig*eig[0][2])
        ten2 = (eig[0][1]/sum_eig)**(-1)*(eig[0][0]/sum_eig*eig[0][2])
        ten3 = (eig[0][2]/sum_eig)**(-1)*(eig[0][0]/sum_eig*eig[0][2])
        a = np.abs(eig[0][0]- eig[0][2])
        b = np.abs(eig[0][1]- eig[0][2])
        c = np.abs(eig[0][1]- eig[0][0])
        ten1 = eig[0][0]+a 
        ten2 = eig[0][1]-c + b
        ten3 = eig[0][2]-a
        eig = (np.array([ten1,ten2,ten3]),eig[1])

    return eig


"""
This function reconstructs a covarinace matirx out of the the object returned
by struc_ten_eig.
The function returns a symmetric matrix, np.array([[xx,xy,xz],
                                                   [xy,yy,yz],
                                                   [xz,yz,zz]])
"""
"""
Function name : revers_eig()
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
def revers_eig(struc_ten):
    rot_mat = struc_ten[1]
    eye = np.eye(3, dtype=np.float64)
    np.fill_diagonal(eye,struc_ten[0])
    
    cov_mat = rot_mat@eye@rot_mat.T    
    return cov_mat


"""
This function takes an unit vector (vector3d), np.array([x,y,z]) and 
calcucates the polar and azimutal angles in respect to 
the Cartesic coordinate system.
Both the vector as well as the angle pair are retuned as np.array([[]]).
This is the format to easly concatenate them to angle_list and vector_list of
any class object in Cone ord NodeEdge.
"""
"""
Function name : convert_vector_angle()
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
def convert_vector_angle(vector3d):
    pol = np.rad2deg(np.arccos(vector3d[2]))
    if pol > 180:
        pol = 360 - pol
    if pol == 0:
        az = 0
    else:
        az = np.rad2deg(np.arctan2(vector3d[1],vector3d[0]))
        if az < 0:
            az = 360 + az
    return np.array([vector3d]), np.array([[pol,az]])


"""
This function is called if a cone decides to branch out. Depending on its 
surroundings (force_struc) and internal branching parameter like the
branching_angle and branching_aperture a branching direction (branching_vector 
and branching_angle) are randomly calculated. The is done via the same Monte
Carlo simultion as for a growth direction.
The paramters are:
    cod_obj : a Growth_Cone Class object
    force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 componetes
                  for the diffusion covarinace
The function returns branch_angle and branch_vector as np.array([[]])
"""
"""
Function name : pdf_convolver_branching()
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
def pdf_convolver_branching(cod_obj, force_drift, force_struc):
    'set new vector'
    diffxx = force_struc[0]
    diffxy = force_struc[1]
    diffxz = force_struc[2]
    diffyy = force_struc[3]
    diffyz = force_struc[4]
    diffzz = force_struc[5]
    struc_ten = np.array([[diffxx,diffxy,diffxz],[diffxy,diffyy,diffyz],[diffxz,diffyz,diffzz]])
    eig = struc_ten_eig(struc_ten, cod_obj.max_eig, cod_obj.min_eig, cod_obj.proxy_reverse_eig)
    cov = revers_eig(eig)
 
    
    coff = np.linalg.solve(eig[1].T,cod_obj.vector_list[-1])
    #print(coff)
    vector_proj2 = coff[0]*eig[1][0]
    vector_proj3 = coff[1]*eig[1][1]
    vector_proj23 = (vector_proj2+vector_proj3)/np.linalg.norm((vector_proj2+vector_proj3))
    vector_proj23, polar_angle = convert_vector_angle(vector_proj23)
    
    montecarlo_iter = cod_obj.montecarlo_iterations
    branching_aperture = 20
    vector_a = fAG.sphere_vector(branching_aperture, polar_angle[0], cov, np.array([0,0,0]), montecarlo_iter, offset =cod_obj.branching_angle[0], points_on_sphere = 4000, plot = False)
    #vector_a = fAG.sphere_vector(branching_aperture, cod_obj.angle_list[-1], cov, np.array([0,0,0]), montecarlo_iter, offset =cod_obj.branching_angle[0], points_on_sphere = 4000, plot = False)
    'new unit vector and its sherical coordinates for growth direction'
    branch_vector_a, branch_angle_a = convert_vector_angle(vector_a)

    return branch_vector_a, branch_angle_a


"""
This function is called if a cone decides to bifurcate. Depending on its 
surroundings (force_struc) and internal bifuration parameter like the
bifucation_angle and bifurcation_aperture two splitting directions 
(branching_vector, branching_angle a and b) are randomly calculated. 
The is done via the same MonteCarlo simultion as for a growth direction.
The paramters are:
    cod_obj : a Growth_Cone Class object
    force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 componetes
                  for the diffusion covarinace
The function returns two a,b branch_angle and branch_vector as np.array([[]])
"""
"""
Function name : pdf_convolver_bifurcation()
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
def pdf_convolver_bifurcation(cod_obj, force_drift, force_struc):
    'set new vector'
    diffxx = force_struc[0]
    diffxy = force_struc[1]
    diffxz = force_struc[2]
    diffyy = force_struc[3]
    diffyz = force_struc[4]
    diffzz = force_struc[5]
    struc_ten = np.array([[diffxx,diffxy,diffxz],[diffxy,diffyy,diffyz],[diffxz,diffyz,diffzz]])
    eig = struc_ten_eig(struc_ten, cod_obj.max_eig, cod_obj.min_eig, cod_obj.proxy_reverse_eig)
    cov = revers_eig(eig)
 
    
    coff = np.linalg.solve(eig[1].T,cod_obj.vector_list[-1])
    #print(coff)
    vector_proj2 = coff[0]*eig[1][0]
    vector_proj3 = coff[1]*eig[1][1]
    vector_proj23 = (vector_proj2+vector_proj3)/np.linalg.norm((vector_proj2+vector_proj3))
    vector_proj23, polar_angle = convert_vector_angle(vector_proj23)
    
    montecarlo_iter = cod_obj.montecarlo_iterations
    bifurcation_aperture = 20
    vector_a = fAG.sphere_vector(bifurcation_aperture, polar_angle[0], cov, np.array([0,0,0]), montecarlo_iter, offset =cod_obj.branching_angle[0], points_on_sphere = 4000, plot = False)
    #vector_a = fAG.sphere_vector(bifurcation_aperture, cod_obj.angle_list[-1], cov, np.array([0,0,0]), montecarlo_iter, offset =cod_obj.bifurcation_angle[0], points_on_sphere = 4000, plot = False)
    vector_b = Quaternions.vec_angle_rot(vector_a,np.array([cod_obj.bifurcation_angle[0]*2]),np.cross(vector_a,cod_obj.vector_list[-1]))
    'new unit vector and its sherical coordinates for growth direction'
    branch_vector_a, branch_angle_a = convert_vector_angle(vector_a)
    vector_b =vector_b+np.random.rand(3)*0.1*np.array([np.random.choice([-1,1]),np.random.choice([-1,1]),np.random.choice([-1,1])])
    vector_b = vector_b/np.linalg.norm(vector_b)
    branch_vector_b, branch_angle_b = convert_vector_angle(vector_b)
                                                           
    return branch_vector_a, branch_angle_a, branch_vector_b, branch_angle_b   

if __name__=='__main__':
    pass
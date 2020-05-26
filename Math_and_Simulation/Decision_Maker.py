#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:02:34 2017

@author: top40ub
"""
import numpy as np


"""the python files containing all needed self writen function"""

import Growth_Classes.Cone_Classes as Cone_Classes

import Math_and_Simulation.Pdf_Convolver as Pdf_Convolver

import Growth_Classes.NodeEdge_Classes as NodeEdge_Classes



"""
Function name : get_deathsilon()
***Description***

The function take as imput a threshold value between 0 and 1 
If a random number d is below the threshold it is set to the Boolean value True
otherwise to the value False

The value d decides if a GrowthCone class object terminates.

***I/O***
Input parameter: 
	a) death_threshold type('array').shape(1,) holds a value between 0 and 1
Output:
	a) d type('boolean') a Boolean to decide if a termination event occurs

Inline output:
Plot output:
Save file:
"""
def get_deathsilon(death_threshold):
    
    np.random.seed()
    d = np.random.random()
    if d < death_threshold:
        d = True          
        d = False #disable termination
    elif d >= death_threshold:
        d = False                                    

    return d


"""
Function name : get_reactsilon()
***Description***

The function take as imput a threshold value between 0 and 1 
If a random number r is below the threshold it is set to the Boolean value True
otherwise to the value False

The value r decides if a terminated GrowthCone class object reactivates.

***I/O***
Input parameter: 
	a) react_threshold type('array').shape(1,) holds a value between 0 and 1
Output:
	a) d type('boolean') a Boolean to decide if a reactivation event occurs

Inline output:
Plot output:
Save file:
"""
def get_reactsilon(react_threshold):
    
    np.random.seed()
    r = np.random.random()
    if r < react_threshold:         
        r = True
        r = False #disable reactivation
    elif r >= react_threshold:
        r = False                                    

    return r


"""
Function name : get_bepsilon()
***Description***

The function take as imput two threshold values each between 0 and 1. 
A set of two random numbers a,b are below their thresholds
(branching_threshold for a and bifurcation_threshold for b) an integer e is set
to 1 for branching or 2 for bifurcation, in case of both event a branching 
occurs, otherwise to the value is 0 for no event.

***I/O***
Input parameter: 
	a) branching_threshold type('array').shape(1,) holds a value 
    between 0 and 1
    b) bifucation_threshold type('array').shape(1,) holds a 
    value between 0 and 1
Output:
	a) e type('int') an integer with value 0,1 or 2 to decide if nothing, a 
    branching event or a bifurcation happen

Inline output:
Plot output:
Save file:
"""
def get_bepsilon(branching_threshold, bifurcation_threshold):
    
    a = np.random.random()
    b = np.random.random()
    np.random.seed()
    
    if a < branching_threshold and b < bifurcation_threshold:
        e = 1
        #e = 0 #disable branching
    elif a < branching_threshold and b >= bifurcation_threshold:
        e = 1
        #e = 0 #disable branching
    elif a >= branching_threshold and b < bifurcation_threshold:
        e = 2
        #e = 0 #disable bifurcation
    else:
        e = 0
    return e


"""
Function name : branching_decision()
***Description***

If a GrowthCone class object decides to branch during one growth step a 
branching is performed.
The function 

***I/O***
Input parameter: 
	a) cod_obj
    b) force_drift
    c) force_struc
    d) view
    e) new_cod
    f) new_nod
    g) new_edd
    h) eddp
    j) nodp
    k) flavourp

Output:
    a) a
    b) new_cod
    c) new_nod
    d) new_edd

Inline output:
Plot output:
Save file:
"""
def branching_decision(cod_obj, force_drift, force_struc, view,
                     new_cod, new_nod, new_edd,
                     eddp, nodp, ed_namesp, no_namesp, co_namesp,
                     flavourp):
    """
    new branching cone
    name for the new branch cone
    """
    cone1 = co_namesp.pop(0)
    a = cod_obj.branching(cone1)
    splitting_mother_cone(cod_obj, new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp)
    nodp[cod_obj.node_list[-2]].edges.append(cod_obj.current_edge)
    

    new_cod[a] = Cone_Classes.BranchingCone(a, *cod_obj.pos_list[-1])
    new_cod[a].node_list.append(cod_obj.node_list[-2])
    branch_vector, branch_angle = Pdf_Convolver.pdf_convolver_branching(cod_obj, force_drift, force_struc)        
    """Starting point parameter"""
    new_cod[a].angle_list = branch_angle
    new_cod[a].vector_list = branch_vector
    new_cod[a].vec_mem = branch_vector[0]
    """
    match internal, memory and growth paramteres of new cone with the old ones
    """
    match_cone_attributes(cod_obj, new_cod[a])
    
    splitting_node_construction(new_cod[a], new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp)
    nodp[cod_obj.node_list[-2]].edges.append(new_cod[a].current_edge)
    
    return a, new_cod, new_nod, new_edd      
    

"""
Function name : bifurcation_decision()
***Description***

--plain text ---

***I/O***
Input parameter: 
	a) cod_obj
    b) force_drift
    c) force_struc
    d) view
    e) new_cod
    f) new_nod
    g) new_edd
    h) eddp
    j) nodp
    k) flavourp
Output:
    a) a
    b) b
    c) new_cod
    d) new_nod
    e) new_edd

Inline output:
Plot output:
Save file:
""" 
def bifurcation_decision(cod_obj, force_drift, force_struc, view,
                         new_cod, new_nod, new_edd,
                         eddp, nodp, ed_namesp, no_namesp, co_namesp,
                         flavourp):
    'new bifurcation cone a, b'
    'name for the new cone a,b'
    cone1 = co_namesp.pop(0)
    cone2 = co_namesp.pop(0)
    a, b = cod_obj.bifurcation(cone1,cone2)
    bifur_vector_a, bifur_angle_a, bifur_vector_b, bifur_angle_b = Pdf_Convolver.pdf_convolver_bifurcation(cod_obj, force_drift, force_struc)
    
    new_cod[a] = Cone_Classes.BifurcationCone(a, *cod_obj.pos_list[-1])
    new_cod[a].node_list.append(cod_obj.node_list[-1])        
    'Starting point parameter'
    new_cod[a].angle_list = bifur_angle_a
    new_cod[a].vector_list = bifur_vector_a
    new_cod[a].vec_mem = bifur_vector_a[0]
    'match internal, memory and growth paramteres of new cone with the old ones' 
    match_cone_attributes(cod_obj, new_cod[a])
    
    splitting_node_construction(new_cod[a], new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp)
    nodp[cod_obj.node_list[-2]].edges.append(new_cod[a].current_edge)
    
    new_cod[b] = Cone_Classes.BifurcationCone(b, *cod_obj.pos_list[-1])
    new_cod[b].node_list.append(cod_obj.node_list[-1])        
    'Starting point parameter'
    new_cod[b].angle_list = bifur_angle_b
    new_cod[b].vector_list = bifur_vector_b
    new_cod[b].vec_mem = bifur_vector_b[0]
    'match internal, memory and growth paramteres of new cone with the old ones' 
    match_cone_attributes(cod_obj, new_cod[b])
    
    splitting_node_construction(new_cod[b],new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp)
    nodp[cod_obj.node_list[-2]].edges.append(new_cod[b].current_edge)
    return a, b, new_cod, new_nod, new_edd


"""
Function name : match_cone_attributes()
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
def match_cone_attributes(cod_obj, new_cod_obj):
    "Internal paramter for new branch"
    new_cod_obj.growth_length = cod_obj.growth_length
    new_cod_obj.aperture = cod_obj.aperture
    new_cod_obj.covariance_revert = cod_obj.covariance_revert
    new_cod_obj.branching_angle = cod_obj.branching_angle
    new_cod_obj.bifurcation_angle = cod_obj.bifurcation_angle

    "Propabilities for branching, bifuraction, termination, reactivation..."
    "... and Monte Carlo iterations"
    new_cod_obj.branchingprob = cod_obj.branchingprob
    new_cod_obj.bifurcationprob = cod_obj.bifurcationprob
    new_cod_obj.deathprob = cod_obj.deathprob
    new_cod_obj.reactivationprob = cod_obj.reactivationprob
    new_cod_obj.montecarlo_iterations = cod_obj.montecarlo_iterations
    "Memory and imformation about nodes, edges and splitting events"
    "and infomation about the compfield/substrate and search field"
    new_cod_obj.searchfield = cod_obj.searchfield
    new_cod_obj.field_dim = cod_obj.field_dim
    new_cod_obj.max_drift = cod_obj.max_drift
    new_cod_obj.min_drift = cod_obj.min_drift
    new_cod_obj.max_eig = cod_obj.max_eig
    new_cod_obj.min_eig = cod_obj.min_eig
    new_cod_obj.last_splitting = 0
    new_cod_obj.memory = cod_obj.memory
    new_cod_obj.proxy_drift = cod_obj.proxy_drift
    new_cod_obj.proxy_tensor = cod_obj.proxy_tensor
    new_cod_obj.proxy_corr = cod_obj.proxy_corr
    new_cod_obj.proxy_reverse_eig = cod_obj.proxy_reverse_eig     


"""
Function name : splitting_node_construction()
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
def splitting_node_construction(cod_obj, new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp):
    'actions on cod_obj'
    'new flavour for branch'
    'new node name for branch'    
    cod_obj.flavour = flavourp.pop(0)
    node = cod_obj.node_construction(no_namesp.pop(0))
    cod_obj.node_list.append(node)
    'new edge name -> the new flavour for branch'
    edge = cod_obj.edge_construction(ed_namesp.pop(0))
    cod_obj.current_edge = [edge]
    
    'new cone node for branch'
    new_nod[node] = NodeEdge_Classes.Node(node,*cod_obj.pos_list[-1])
    new_nod[node].constructor_cone.append(cod_obj.name)
    new_nod[node].edges.append(edge)
    
    'new edge for branch'
    new_edd[edge] = NodeEdge_Classes.Edge(edge,*cod_obj.pos_list[-1],*cod_obj.pos_list[-1])
    new_edd[edge].constructor_cone.append(cod_obj.name)
    new_edd[edge].nodes.append(cod_obj.node_list[-2])
    new_edd[edge].nodes.append(node)
    new_edd[edge].pos_list = np.array([cod_obj.pos_list[-1]])
    new_edd[edge].flavour = cod_obj.flavour


"""
Function name : splitting_mother_cone()
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
def splitting_mother_cone(cod_obj, new_nod, new_edd, eddp, nodp, ed_namesp, no_namesp, co_namesp, flavourp):
    'actions on cod_obj'
    'new node name for branch'    
    node = cod_obj.node_construction(no_namesp.pop(0))
    cod_obj.node_list.append(node)
    'new edge name -> the new flavour for branch'
    edge = cod_obj.edge_construction(ed_namesp.pop(0))
    cod_obj.current_edge = [edge]
    cod_obj.last_splitting = 0
        
    'new cone node for branch'
    new_nod[node] = NodeEdge_Classes.Node(node,*cod_obj.pos_list[-1])
    new_nod[node].constructor_cone.append(cod_obj.name)
    new_nod[node].edges.append(edge)
    #new_nod[node].edges.append(cod_obj.current_edge)
    'new edge for branch'
    new_edd[edge] = NodeEdge_Classes.Edge(edge,*cod_obj.pos_list[-1],*cod_obj.pos_list[-1])
    new_edd[edge].constructor_cone.append(cod_obj.name)
    new_edd[edge].nodes.append(cod_obj.node_list[-2])
    new_edd[edge].nodes.append(node)
    new_edd[edge].pos_list = np.array([cod_obj.pos_list[-1]])
    new_edd[edge].flavour = cod_obj.flavour
    


if __name__=='__main__':
    pass     
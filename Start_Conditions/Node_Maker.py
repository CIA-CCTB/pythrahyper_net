#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 16:38:25 2017

@author: top40ub
"""

import numpy as np
from pathlib import Path
import Growth_Classes.Cone_Classes as Cone_Classes 
import Growth_Classes.NodeEdge_Classes as NodeEdge_Classes


""" This function constructs the first seeds/nodes of the growth cones """
"""
Function name : cone_construction()
***Description***

The function constructs the initial GrowthCone class objects. 
Each cone is initialised by its name 'cn' for each of the n starting positions
and start orientations of the imput array start_pos and st_angle.
All needed parametes for dem paramter dictionaries internal, growth, multilayer
are transfered to each GrowthCone class object. Each GrowthCone class object
collects its own flavour from the flavor list.

A dictionary cod with key = 'cn' with for each for the n positions
and value = GrowthCone class object for each position is returned

***I/O***
Input parameter: 
	a) start_pos type(nd.array).shape(n,3) array with all 3D starting postions 
	of the cone 
	b) st_angle type(nd.array).shape(n,2) array with all starting orientations 
	of the cone (pol,az)
	c) internal type('dict') dictionary with all internal GrowthCone class and
	NodeEdge class object parameter
	d) growth type('dict') dictionary with all parameter for the growth simulation
	e) multilayer type('dict') dictionary with all parameter of the multilayer
	computation field
	f) flavour type('list') the flavour list to distribute a flavour to each 
	GrowthCone class object
Output:
	a) cod type('dict') dictionary with all n starting Cone Class objects

Inline output:
Plot output:
Save file:
"""
def cone_construction(start_pos, st_angle, internal, growth, multilayer, flavour, cone_names):

    cone_start_list = {}
    st_angle_list = {}
    if len(start_pos) != len(st_angle):
        raise Exception("Number of start positions don't match with number of start angles")
    'define grothcone dictionary [cod], growth aperture [aperture]'
    cod = {}
    for i in range(len(start_pos)):
        if len(start_pos[i]) != 3:
            raise Exception("Three dimensional vector for each start postion is needed (x,y,z)")
        if len(st_angle[i]) != 2:
            raise Exception("Two angles (pol,az) are needed for each the start orientation")
        j = cone_names.pop(0)
        cone_start_list['c' + str(j)] =start_pos[int(j-1)]
        st_angle_list['c' + str(j)] = st_angle[int(j-1)]
        del i

    for key in cone_start_list.keys():
        cod[key] = Cone_Classes.GrowthCone(key, cone_start_list[key][0], cone_start_list[key][1], cone_start_list[key][2])
        
        "Starting point parameter"
        cod[key].angle_list = np.array([st_angle_list[key]])
        pol = cod[key].angle_list[0][0]
        az  = cod[key].angle_list[0][1]
        cod[key].vector_list = np.array([[np.sin(np.deg2rad(pol)) * np.cos(np.deg2rad(az)), np.sin(np.deg2rad(pol)) * np.sin(np.deg2rad(az)), np.cos(np.deg2rad(pol))]])
        cod[key].vec_mem = cod[key].vector_list[0]
        "Internal paramter"
        cod[key].growth_length = internal['stepsize']
        cod[key].aperture = internal['aperture']
        cod[key].covariance_revert = growth['covariance_revert']
        cod[key].branching_angle = internal['branching_angle']
        cod[key].bifurcation_angle = internal['bifurcation_angle']
        cod[key].memory = internal['memory']

        "Propabilities for branching, bifuraction, termination, reactivation..."
        "... and Monte Carlo iterations"
        cod[key].branchingprob = internal['branchingprob']
        cod[key].bifurcationprob = internal['bifurcationprob']
        cod[key].deathprob = internal['deathprob']
        cod[key].reactivationprob = internal['reactivationprob']
        cod[key].montecarlo_iterations = growth['montecarlo_iterations']
        
        "Memory and imformation about nodes, edges and splitting events"
        "and infomation about the compfield/substrate and search field"
        cod[key].searchfield = internal['searchfield'].reshape(3,3,3)
        cod[key].field_dim = multilayer['dim']
        cod[key].max_drift = multilayer['max_drift']
        cod[key].min_drift = multilayer['min_drift']
        cod[key].max_eig = multilayer['max_eigenvalue']
        cod[key].min_eig = multilayer['min_eigenvalue']
        cod[key].frame_seq.append(0)
        cod[key].flavour = flavour.pop(0)
        """Proxy_PDF_Parameter"""
        cod[key].proxy_drift = growth['Proxy_drift']
        cod[key].proxy_tensor = growth['Proxy_tensor']
        cod[key].proxy_corr = growth['Proxy_corr']
        cod[key].proxy_reverse_eig = growth['Proxy_reverse_eig']           
        del key
    return cod

"""
Function name : node_construction()
***Description***

The function creates for each unique element in the common_node
array a NodeEdge class object (node) with the starting point for the node at
the position in the start_pos array for the first unique element in 
common_node.
The node name is 'n0',...,'nn' where n indicates the unique elements.

A dictionary cod with key = 'n0',...'nn' with for each for the n uniques
and value = NodeEdge class object (node) for each unique position is returned

***I/O***
Input parameter: 
	a) start_pos type(np.'ndarray').shape(n,3) array with all 3D starting postions 
	of the cone
	b) common_node type('np.ndarray').shape(n,1) array with the name of the starting node
	for each GrowthCone class object
Output:
	a) nod type('dict') dictionary with all n starting NodeEdge class objects (node)

Inline output:
Plot output:
Save file:
"""
def node_construction(start_pos, common_node, node_names):
    node_start_list ={}
    nod = {}
    uni, ind_uni = np.unique(common_node, return_index = True)
    for i, ind in zip(uni, ind_uni):
        j = node_names.pop(0)
        node_start_list['n' + str(int(j))] = start_pos[int(ind)]
    for key in node_start_list.keys():
        nod[key] = NodeEdge_Classes.Node(key, node_start_list[key][0], node_start_list[key][1], node_start_list[key][2])    
    return nod


"""
Function name : cod_nod_edd_match()
***Description***

The function constructs the first edge for each of of the n GrowthCone class objects.
Wit the information in cod, nod and the common_node array all parameter and references
between the three class objects are matched and updated.

Each cone points towards two nodes, its starting node and one at its tip that moves 
alongside.
Each cone has starts with on edge connection those nodes. During the growth the edge is
elongated. These edges inherited the flavour of the cone.
Each node lists all edges connecting to it.

Three dictionaries cod, nod, edd for the GrowthCone NodeEdge (node), NodeEdge (edge) class
objects are returned

***I/O***
Input parameter: 
	a) cod type('dict') dictionary with all n starting Cone Class objects
	b) nod type('dict') dictionary with all n starting NodeEdge class objects (node)
	c) common_node type('np.ndarray').shape(n,1) array with the name of the starting node
	for each GrowthCone class object
Output:
	a) cod type('dict') dictionary with all n starting Cone Class objects updated
	b) nod type('dict') dictionary with all n starting NodeEdge class objects (node) updated
	c) edd type('dict') dictionary with all starting edges

Inline output:
Plot output:
Save file:
"""
def cod_nod_edd_match(cod, nod, common_node, node_names, edge_names):
    edd = {}
    new_nod = {}
    for keytup in zip(cod,common_node):
        nod['n' + str(int(keytup[1]))].constructor_cone.append(keytup[0])
        cod[keytup[0]].node_list.append('n' + str(int(keytup[1])))
        
        node2 = cod[keytup[0]].node_construction(node_names.pop(0))
        new_nod[node2] = NodeEdge_Classes.Node(node2,*cod[keytup[0]].pos_list[0])
        
        cod[keytup[0]].node_list.append(node2)
        edge =cod[keytup[0]].edge_construction(edge_names.pop(0))
        cod[keytup[0]].current_edge = [edge]
        nod['n' + str(int(keytup[1]))].edges.append(edge)
        new_nod[node2].edges.append(edge)
        edd[edge] = NodeEdge_Classes.Edge(edge,*cod[keytup[0]].pos_list[0],*cod[keytup[0]].pos_list[0])
        edd[edge].flavour=cod[keytup[0]].flavour
        edd[edge].constructor_cone.append(keytup[0])
        edd[edge].nodes.append('n' + str(int(keytup[1])))
        edd[edge].nodes.append(node2)
        edd[edge].pos_list=cod[keytup[0]].pos_list
    nod.update(new_nod)
    return cod, nod, edd



if __name__=='__main__':
    pass
    

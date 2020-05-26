#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:40:57 2020

@author: top40ub
"""

import numpy as np
import csv
from pathlib import Path
import os

"""
Function name : import_internal_par()
***Description***

The function loads and reads the .csv-file withthe internal parameters
for the Cone and NodeEdge class objects. The imput variable of type('str')
contains then name and path to the parameter file.
A dictionary with key, value pairs for the paramters is returned

***I/O***
Input parameter: 
	a) path_internal type('str') path to the .csv-file containing all
	internal parameter for the Cone class and NodeEdge class objects
Output:
	a) internalpardict type('dict') key, value pairs for all internal 
	parameters

Inline output:
Plot output:
Save file:
"""
def import_internal_par(path_internal):
    internalpardict = {}
    with open(path_internal, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=',')
        L = list(filereader)
        while [] in L:
            L.pop(L.index([]))
        for row in L:
            if len(row) >= 2:
                internalpardict[row[0]] = np.array(row[1:], dtype=float)

    return internalpardict


"""
Function name : import_growth_par()
***Description***

The function loads and reads the .csv-file with the growth parameters
for the growth simulation. The imput variable of type('str')
contains then name and path to the parameter file.
A dictionary with key, value pairs for the paramters is returned.

***I/O***
Input parameter: 
	a) path_growth type('str') path to the .csv-file containing all
	growth parameter for the growth simulation
Output:
	a) growthpardict type('dict') key, value pairs for all growth 
	parameters

Inline output:
Plot output:
Save file:
"""
def import_growth_par(path_growth):
    growthpardict = {}
    with open(path_growth, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=',')
        L = list(filereader)
        while [] in L:
            L.pop(L.index([]))
        for row in L:
            if len(row) >= 2:
                if row[1:][0] == 'True':
                    growthpardict[row[0]] = True
                elif row[1:][0] == 'False':
                    growthpardict[row[0]] = False
                else:
                    growthpardict[row[0]] = np.array(row[1:], dtype=float)

    return growthpardict


"""
Function name : import_multilayer_par()
***Description***

The function loads and reads the .csv-file with the multilayer parameters
for the growth simulation. The imput variable of type('str')
contains then name and path to the parameter file.
A dictionary with key, value pairs for the paramters is returned.

***I/O***
Input parameter: 
	a) path_multilayer type('str') path to the .csv-file containing all
	parameter of the multilayer computation field
Output:
	a) multilayerdict type('dict') key, value pairs for all multilayer 
	parameters

Inline output:
Plot output:
Save file:
"""
def import_multilayer_par(path_multilayer):
    multilayerdict = {}
    with open(path_multilayer, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=',')
        L = list(filereader)
        while [] in L:
            L.pop(L.index([]))
        for row in L[1:6]:
            multilayerdict[row[0]] = np.array(row[1:], dtype=float)
        for row in L[6:]:
            multilayerdict[row[0]] = row[1:]
    return multilayerdict


"""
Function name : import_startpos_par()
***Description***

The function loads and reads the .csv-file with the start position
parameter for the growth simulation. The imput variable of type('str')
contains then name and path to the parameter file.
Three ndarrays with shape (n,3) for all start position (n,2) start orientations
and (n,1) starting node are returned.

***I/O***
Input parameter: 
	a) path_startpar type('str') path to the .csv-file containing the position
	and orientation of all starting Cone class objects and all NodeEdge class
	objects
Output:
	a) np.array(pos_list,dtype=float).shape(n,3) start postion of each cone
	b) np.array(angle_list, dtype=float).shape(n,2) start orientation
	c) np.array(common_node, dtype=float).shape(n,1) starting node

Inline output:
Plot output:
Save file:
"""
def import_startpos_par(path_startpar):
    pos_list = []
    angle_list = []
    common_node = []
    with open(path_startpar, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=',')
        L = list(filereader)
        while [] in L:
            L.pop(L.index([]))
        for i in L[2:]:
            pos_list += [i[:3]]
            angle_list += [i[3:5]]
            common_node += [i[5:]] 
    return np.array(pos_list,dtype=float), np.array(angle_list, dtype=float), np.array(common_node, dtype=float)

if __name__=='__main__':
    home = str(Path(os.getcwd()).parent)
    path_internal = home + '/Parameter_Files/internal_parameters.csv'
    path_growth = home + '/Parameter_Files/growth_parameters.csv'
    path_multilayer = home + '/Parameter_Files/multilayer_dir_parameters.csv'
    path_startpar = home + '/Parameter_Files/starting_positions.csv'
    d1 = import_internal_par(path_internal)
    d2 = import_growth_par(path_growth)
    d3 = import_multilayer_par(path_multilayer)
    d40,d41,d42 = import_startpos_par(path_startpar)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:17:37 2020

@author: top40ub
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from mpl_toolkits.mplot3d import Axes3D


def cones_within_region(cod, box):
    L = []
    for key in cod.keys():
        r = np.sqrt(np.sum((cod[key].pos_list-125)**2, axis= 1))
        x,y,z = cod[key].pos_list.T
        if r[r>=99].size != 0:
            L.append(key)
    print(L)
    return L           

def edge_within_region(edd, box):
    L = []
    for key in edd.keys():
        for i in edd[key].pos_list:
            x_pos, y_pos, z_pos = i
            if (x_pos <= 30 or x_pos >= 220
                or y_pos <= 30 or y_pos >= 220
                or z_pos <= 30 or z_pos >= 220):
                L.append(key)
 return L

def traceback_cod_splitting(key, cod, nod, edd):
    current_edge = cod[key].current_edge
    pass


"""
Function name : create_edd_length()
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
def create_edd_length(edd, plot = True):
    edge_len_dis = []
    for key in edd.keys():
        n_n_dis = np.linalg.norm(edd[key].pos_list[0]-edd[key].pos_list[-1])
        curve_len = np.array([0.])
        for i in edd[key].pos_list:
            curve_len += np.linalg.norm(i) 
        steps = edd[key].pos_list.shape[0]
        edge_len_dis += [n_n_dis]
    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(edge_len_dis, 40, density=1, facecolor='green', alpha=0.75)
        plt.show()
    return n_n_dis, curve_len, steps 


if __name__=='__main__':
    n_n_dis, curve_len, steps = create_edd_length(edd)
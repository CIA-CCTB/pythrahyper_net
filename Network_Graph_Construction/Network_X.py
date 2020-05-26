#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:30:48 2018

@author: top40ub
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from mpl_toolkits.mplot3d import Axes3D



"""
Function name : create_graph_from_edd_nod()
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
def create_graph_from_edd_nod(edd, nod, plot_net = True, plot_net_circ = True, plot_net_spring = True):
 
    G = nx.Graph()
    
    nodes = [key for key in nod.keys()]
    
    G.add_nodes_from(nodes)
    
    edges = [(edd[key].nodes[0],edd[key].nodes[1]) for key in edd.keys()]
    
    G.add_edges_from(edges)
    subG = list(nx.connected_components(G))
    
    if plot_net == True:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        nx.draw_networkx(G,ax = ax1,with_labels = False,node_size = 10, node_color =  'b', node_shape = 'o',width = 0.5)
    if plot_net_circ == True:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        nx.draw_circular(G,ax = ax2,with_labels =False,node_size = 10, node_color =  'b', node_shape = 'o',width = 0.5)
    if plot_net_spring == True:
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        nx.draw_spring(G,ax = ax3,with_labels = False,node_size = 10, node_color =  'b', node_shape = 'o',width = 0.5)    
    
    plt.draw()
    
    return G, subG


   
"""
Function name : voronoi_tesselation()
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
def voronoi_tesselation(cod, dim, plot_ridge_line = True):
        points = np.array([[0,0,0]])
        for key in cod:
            points = np.concatenate((points,[cod[key].pos_list[-1]]),axis=0)
    
        vor = Voronoi(points)
        fig = plt.figure()
        ax3d = fig.add_subplot(111, projection='3d')
        if plot_ridge_line == True:
            for i in vor.ridge_vertices:
                if -1 not in i:
                    if np.amin(vor.vertices[i]) < 0 or np.amax(vor.vertices[i]) > np.amax(dim):
                        pass
                    else:
                        ax3d.plot([vor.vertices[i[k]][0] for k in range(len(i))],
                                     [vor.vertices[i[k]][1] for k in range(len(i))],
                                     [vor.vertices[i[k]][2] for k in range(len(i))], c = 'r',linewidth=0.5)
        for key in cod.keys():
            ax3d.scatter(*cod[key].pos_list[-1],c='b')
        ax3d.set_xlim([0, dim[0]])
        ax3d.set_ylim([0, dim[0]])
        ax3d.set_zlim([0, dim[0]])
        plt.show()
        
        return vor, fig, ax3d


"""
Function name : cone_network_cluster()
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
def cone_network_cluster(cod,radius, plot_network = True):

    G = nx.Graph()
    for key in cod:
        G.add_node(key, pos = cod[key].pos_list[-1])
        
    for key1 in G.nodes:
        for key2 in G.nodes:
            if key1 > key2 and np.linalg.norm(abs(np.array(G.nodes[key1]['pos'])-np.array(G.nodes[key2]['pos']))) < radius:
                G.add_edge(key1,key2)
    
    subG = list(nx.connected_components(G))
    if plot_network == True:
        fig_cluster, ax_cluster = plt.subplots()
        
        options = {'node_color': 'red',
                   'node_size': 4,
                   'edge_color': 'green',
                   'width': 2} 
        nx.draw(G,with_labels =False,**options)       
        plt.show() 

    return G, subG


"""
Function name : voronoi_network_cluster()
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
def voronoi_network_cluster(cod, dim, radius, plot_network = True, plot_ridge_line = True):
    
    G_vor, subG_vor = cone_network_cluster(cod,radius, plot_network= plot_network)
    vor, fig, ax3d = voronoi_tesselation(cod, dim, plot_ridge_line= plot_ridge_line)
    
    
    for i in list(G_vor.edges):
        ax3d.plot([cod[i[0]].pos_list[-1][0],cod[i[1]].pos_list[-1][0]],
                  [cod[i[0]].pos_list[-1][1],cod[i[1]].pos_list[-1][1]],
                  [cod[i[0]].pos_list[-1][2],cod[i[1]].pos_list[-1][2]],
                  c = 'g', alpha = 0.4,linewidth=1)   


    for key in cod.keys():
        
        u = np.linspace(0, 2 * np.pi, 10)
        v = np.linspace(0, np.pi, 10)
        x = radius * np.outer(np.cos(u), np.sin(v)) + cod[key].pos_list[-1][0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + cod[key].pos_list[-1][1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cod[key].pos_list[-1][2]
        ax3d.plot_surface(x, y, z, color='b', alpha= 0.1)
    plt.show()

    """
    vor.regions = [cut for cut in vor.regions if cut != []]
    for c in range(len(subG)):
        for p in subG[c].nodes:
            ind = vor.points.tolist().index(cod[p].tolist())+1
            ind_r = vor.point_region.tolist().index(ind)            
            region = vor.regions[ind_r]
            if not -1 in region:
                polygon = [vor.vertices[j] for j in region]           
                plt.fill(*zip(*polygon), color=(c / len(subG),c / len(subG), c / len(subG)))
    """
    return vor, G_vor, subG_vor

if __name__=='__main__':
    G, subG = create_graph_from_edd_nod(edd,nod)
    n_n_dis, curve_len, steps = create_edd_length(edd)
    vor, G_vor, subG_vor = voronoi_network_cluster(cod, dim, 10)
    vor, fig, ax3d = voronoi_tesselation(cod, dim, plot_ridge_line = True)
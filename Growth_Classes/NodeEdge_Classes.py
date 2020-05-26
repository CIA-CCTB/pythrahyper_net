#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 14:26:30 2018

@author: top40ub
"""

import numpy as np


"""
Class name : NodeEdgeClass()
***Description***

The main NodeEdgeClass class, all further NOde or Edge class objects are subclasses
of this class.
IMPORTANT the main class has no __init__ function.
This class also helds the __del__ function to remove the class object.

***I/O***
Input parameter: 
Output:

Inline output:
Plot output:
Save file:
"""  
class NodeEdgeClass:
    node_counter = 0
    edge_counter = 0

    
    def __del__(self):
        NodeEdgeClass.node_counter -= 1
        NodeEdgeClass.edge_counter -= 1
        
"""
Class name : Node()
***Description***

Subclass to the NodeEdgeClass class it has a __init__ function
that initializes a node at position x,y,z

***I/O***
Input parameter: 
Output:

Inline output:
Plot output:
Save file:
"""         
class Node(NodeEdgeClass):
    '__init__ initializes a certain node at pos x,y,z'

    def __init__(self, name,  x_pos = 0, y_pos = 0, z_pos = 0):
        NodeEdgeClass.node_counter += 1

        self.name = name

        'Starting point of the node'
        self.x_pos = x_pos

        self.y_pos = y_pos

        self.z_pos = z_pos
        
        'pos_list is atribute with an array containing the position of the node in 3d space'
        self.pos_list = np.array([[self.x_pos, self.y_pos, self.z_pos]])
        
        'construction_cone contains the name/keyword of the cone building this node'
        'if the node is a init node or it occurs during branching, bifurcation or terminating the mothercone is the constructor'
        'if the node occurs during merging the cone fusioning with the edge is the constructor'
        self.constructor_cone = []
        
        'edges is a list containing the name//key of all edges connected to the node'
        self.edges = []


"""
Class name : Edge()
***Description***


Subclass to the NodeEdgeClass class it has a __init__ function
that initializes a node at initializes a certain edge between 
two nodes at pos x1,y1,z1 and x2,y2,z2

***I/O***
Input parameter: 
Output:


Inline output:
Plot output:
Save file:
"""  
class Edge(NodeEdgeClass):
    '__init__ initializes a certain edge between two nodes at pos x1,y1,z1 and x2,y2,z2'
    
    def __init__(self, name,  x_pos1 = 0, y_pos1 = 0, z_pos1 = 0, x_pos2 = 0, y_pos2 = 0, z_pos2 = 0 ):
        NodeEdgeClass.edge_counter += 1

        self.name = name
        self.flavour = 0
        'Position of node2'
        self.x_pos1 = x_pos1

        self.y_pos1 = y_pos1

        self.z_pos1 = z_pos1
        
        
        'Position of node2'
        self.x_pos2 = x_pos2
        
        self.y_pos2 = y_pos2

        self.z_pos2 = z_pos2        
        
        'pos_list is atribute with an array containing the position of the nodes in 3d space'
        self.pos_list = np.array([[self.x_pos1, self.y_pos1, self.z_pos1],[self.x_pos2, self.y_pos2, self.z_pos2]])
        
        'construction_cone contains the name/keyword of the cone building this node'
        'if the node is a init node or it occurs during branching, bifurcation or terminating the mothercone is the constructor'
        'if the node occurs during merging the cone fusioning with the edge is the constructor'
        self.constructor_cone = []
        
        'edges is a list containing the name//key of the two nodes connected bye the edge'
        self.nodes = []
        self.bbox()


    """
    Function name : bbox()
    ***Description***
    
    Creates all 8 coner points of a box containing a NodeEdge class object (edge).
    The boundingbox is saved as an class attribute
    
    ***I/O***
    Input parameter: 
    Output:
    
    Inline output:
    Plot output:
    Save file:
    """      
    def bbox(self):
        points= self.pos_list
        a = np.zeros((3,2))
        a[0,0] = np.min(points[:,0])
        a[0,1] = np.max(points[:,0])
        a[1,0] = np.min(points[:,1])
        a[1,1] = np.max(points[:,1])
        a[2,0] = np.min(points[:,2])
        a[2,1] = np.max(points[:,2])
        self.boundingbox = a

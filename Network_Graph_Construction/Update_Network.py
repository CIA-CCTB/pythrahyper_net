#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 00:22:10 2018

@author: top40ub
"""
import numpy as np
import Growth_Classes.NodeEdge_Classes as NodeEdge_Classes


"""
Function name : update_network_edd()
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
def update_network_edd(cod, edd, nod, edge_names, node_names, cone_names):
    merge_tupel_list = []
    new_edd ={}
    key_l = []
    
    for key in cod.keys():
        if len(cod[key].sleep_merge_event[0]) == 4:
            merge_tupel_list.append(cod[key].sleep_merge_event[-1])

    merge_tupel_list =sorted(merge_tupel_list,key=getkey0)

    for key in edd.keys():
        d = []
        node = edd[key].nodes[0]
        pos_start = 0
        cones_merged_with_edge = 0
        
        for tupel in merge_tupel_list:
            if tupel[0] == key:
                d.append(tupel)
                cones_merged_with_edge += 1
                
        if d != []:
            d= sorted(d,key = getkey1)
            for j, c in zip(d,range(cones_merged_with_edge)):
                'new edge that start at pos_start/j[4]'
                new_edge_name = 'e'+ str(edge_names.pop(0))
                pos_end = j[1]
                new_edd[new_edge_name] = NodeEdge_Classes.Edge(new_edge_name, *edd[key].pos_list[0], *cod[j[2]].pos_list[-1])
                new_edd[new_edge_name].constructor_cone.append(j[2])
                new_edd[new_edge_name].nodes.append(node)
                new_edd[new_edge_name].nodes.append(j[3])
                new_edd[new_edge_name].pos_list=edd[key].pos_list[pos_start:pos_end+1]
                new_edd[new_edge_name].flavour = edd[key].flavour
                new_edd[new_edge_name].x_pos2 = edd[key].pos_list[j[1]][0] 
                new_edd[new_edge_name].y_pos2 = edd[key].pos_list[j[1]][1] 
                new_edd[new_edge_name].z_pos2 = edd[key].pos_list[j[1]][2]
                
                pos_start = pos_end
                node = j[3]
                
            cone = edd[key].constructor_cone[-1]
            l_edge = edd[key].name
            #print(pos_end, edd[key].pos_list.shape,l_edge)
            new_edd[l_edge] = NodeEdge_Classes.Edge(l_edge,*edd[key].pos_list[pos_end],*edd[key].pos_list[-1])
            new_edd[l_edge].constructor_cone.append(cone)
            new_edd[l_edge].nodes.append(node)
            new_edd[l_edge].nodes.append(edd[key].nodes[-1])
            cod[cone].current_edge = [l_edge]
            new_edd[l_edge].flavour = edd[key].flavour
            key_l.append(key)
            
        else:
            pass
    edd.update(new_edd)  


"""
Function name : update_network_nod()
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
def update_network_nod(cod, edd, nod, edge_names, node_names, cone_names):
    for key in edd.keys():
        for node in edd[key].nodes:
            if key not in nod[node].edges:
                nod[node].edges.append(key)


"""
Function name : view_flavour_field()
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
'Functions called in sorted operation to sort after item positions'
def getkey0(item):
    return item[0]
def getkey1(item):
    return item[1]
def getkey2(item):
    return item[2]
def getkey3(item):
    return item[3]


if __name__=='__main__':
    pass
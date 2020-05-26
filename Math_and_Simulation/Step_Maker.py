#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 16:55:52 2017

@author: top40ub
"""

#import numpy as np


import numpy as np
from random import shuffle
import multiprocessing


"""
The python file containing all needed self writen functions to deciced the
actions a growth cone can do in one growth step
"""

"""
This python file is a rng to decide if branching, bifurcation or
termination happens
"""
import Math_and_Simulation.Decision_Maker as Decision_Maker
"""
Python file which calls the 3dpdf_maker function but also converts the local
structure tensor into pdfs´ for navigation
"""
import Math_and_Simulation.Pdf_Convolver as Pdf_Convolver


"""
Function name : growing_function()
***Description***

The function is a decision tree that checks the state (sleep_merge_event)
of the growth cone (cod_obj) and its surroundings force_drift, force_struc,
view, eddp, nodp, flavoup, fieldp. Based on this and its internal decision 
properties it
either changes the state cod_obj.sleep_merge_event from [('growing')] 
                                                     to [('sleeping')],
                                                     to [('out of bounds')],
                                                     to [merge_event]
or it stays [('growing')] or reactivates from [('sleeping')] to [('growing')].

If it stays [('growing')] it decides to simple grow (simple_growth)
                                     or to branch (branching_growth)
                                     or to bifurcate (bifurcation_growth).

The growth functions shift all GrowthCone class objects to a new position:
    simple_growth: shifts cod_obj
    branching_growth: shifts cod_obj and new_cod[a]
    bifurction_growth: shifts new_cod[a], new_cod[b] and sets state of cod_obj
                       to [('bifurcation')].
    
    The manager objects are updated:
    eddp and nodp are updated with new_edd and new_nod to match all NodeEdge
    class objects.
    Each NodeEdge Class object (edge) places its flavour reference to its 
    new position in the flavourfield (fieldp) 

***I/O***
Input:
    a) cod_obj : the GrowthCone Class object
    b) force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    c) force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 structure
                  componentes of the substrate
    d) view : a np.array().shape(3,3,3) with information about the next NodeEdge 
           class neighbours
    e) new_cod, new_edd, new_nod : empty dictionaries to intermediate collect 
                                newly constructed GrowthCone class objects and 
                                NodeEdge class objects
    f) eddp : Manager dictionary for NoedEdge class objectes (edge) 
    g) nodp : Manager dictionary for NoedEdge class objectes (node) 
    h) flavourp : Manager shared memory object (list) with the free to .pop(0) 
               flavour clolours there is a max flavour set in beginning
    i) fieldp : Manager shared memory object (array) the flavourfield filled 
             elementwise with references type(int) to the flavour of
             NodeEdge class objects. Each edge position is flavour marked in
             the array. The array is a flatten version 
             (1D reshape(dimx,dimy,dimz)) of the field.
    j) step : the number of the current iteration step
Output:
    a) cob_obj : the GrowthCone Class object in new state/position
    b) new_cod : dictionary containing eventual key = 'name' and GrowthCone class
               object of all newly constructed GrowthCone class objects
    c) new_nod : dictionary containing eventual key = 'name' and NodeEdge class
              object of all newly constructed NodeEdge class objects (node)
    d) new_edd : dictionary containing eventual key = 'name' and NodeEdge class
              object of all newly constructed NodeEdge class objects (edge)

Inline output:
Plot output:
Save file:
"""
def growing_function(cod_obj, force_drift, force_struc, view, 
                     new_cod, new_nod, new_edd,
                     eddp, nodp,
                     ed_namesp, no_namesp, co_namesp, flavourp, fieldp, frame):

    "Recall last position coordinates"
    x_pos, y_pos, z_pos = cod_obj.pos_list[-1]

    "Create decision signal varibales d, r, e, l"
    "Boolian termintation signal"
    "Boolian reactivation signal"
    "Branching and bifuraction probapility 0, 1 or 2"
    "Last branching event"
    d = Decision_Maker.get_deathsilon(cod_obj.terminationprob)
    r = Decision_Maker.get_reactsilon(cod_obj.reactivationprob)
    e = Decision_Maker.get_bepsilon(cod_obj.branchingprob, cod_obj.bifurcationprob)
    l = cod_obj.last_splitting
    
    "Cone information about the field dimensions"
    dim = cod_obj.field_dim
    "Ceck if cone edge is long enough to merge"
    if l <= 3:
        merge_event = False
    else:
        'check if cone can merge/collision detection'
        merge_event = merge_edges(cod_obj, view, eddp)
        #merge_event = False #disable merging
    
    "Check if growth cone is inactivate by checking..."
    "...the sleep_merge attribut"
    "Check if growth cone is moved outside..."
    "...the compfield if this is TRUE it is inactivated"
    "Check if growth cone randomly terminates"

    
    if (x_pos <= 2 or x_pos >= dim[0]-2
        or y_pos <= 2 or y_pos >= dim[1]-2
        or z_pos <= 2 or z_pos >= dim[2]-2):

        cod_obj.sleep_merge_event = [('out_of_bounds')]
        cod_obj.frame_seq.append(frame)
        pass

    elif d == True:
        "decide if cone termiates"
        cod_obj.sleep_merge_event = [('sleeping')]
        cod_obj.frame_seq.append(frame)
        pass

    elif merge_event != False:
        'merge cone with edge if collions happend'
        cod_obj.sleep_merge_event = [merge_event]
        merge_cone_with_edge(cod_obj, merge_event, new_edd, new_nod, eddp, nodp, fieldp)
        cod_obj.frame_seq.append(frame)
        pass
 
    elif cod_obj.sleep_merge_event == [('sleeping')] and r == True:
        "Check if a terminated growth cone reactivates"
        cod_obj.sleep_merge_event = [('growing')]
        cod_obj.frame_seq.append(frame)
        if e == 0:
            simple_growth(cod_obj, force_drift, force_struc, view,
                          new_cod, new_nod, new_edd,
                          eddp, nodp, ed_namesp, no_namesp, co_namesp,
                          flavourp, fieldp, frame)
        
        elif e == 1:
            branching_growth(cod_obj, force_drift, force_struc, view,
                             new_cod, new_nod, new_edd,
                             eddp, nodp, ed_namesp, no_namesp, co_namesp,
                             flavourp, fieldp, frame)

        elif e == 2:
            cod_obj.sleep_merge_event = [('bifurcation')]
            cod_obj.frame_seq.append(frame)
            bifurcation_growth(cod_obj, force_drift, force_struc, view,
                               new_cod, new_nod, new_edd,
                               eddp, nodp, ed_namesp, no_namesp, co_namesp,
                               flavourp, fieldp, frame)
            
    else:
        if e == 0:
            simple_growth(cod_obj, force_drift, force_struc, view,
                          new_cod, new_nod, new_edd,
                          eddp, nodp, ed_namesp, no_namesp, co_namesp,
                          flavourp, fieldp, frame)
            
        elif e == 1:
            branching_growth(cod_obj, force_drift, force_struc, view,
                             new_cod, new_nod, new_edd,
                             eddp, nodp, ed_namesp, no_namesp, co_namesp,
                             flavourp, fieldp, frame)

        elif e == 2:
            cod_obj.sleep_merge_event = [('bifurcation')]
            cod_obj.frame_seq.append(frame)
            bifurcation_growth(cod_obj, force_drift, force_struc, view,
                               new_cod, new_nod, new_edd,
                               eddp, nodp, ed_namesp, no_namesp, co_namesp,
                               flavourp, fieldp, frame)

    return cod_obj, new_cod, new_nod, new_edd


"""
Function name : simple_growth()
***Description***

--plain text ---

***I/O***
Input parameter: 
    a) cod_obj : the GrowthCone Class object
    b) force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    c) force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 structure
                  componentes of the substrate
    d) view : a np.array().shape(3,3,3) with information about the next NodeEdge 
           class neighbours
    e) new_cod, new_edd, new_nod : empty dictionaries to intermediate collect 
                                newly constructed GrowthCone class objects and 
                                NodeEdge class objects
    f) eddp : Manager dictionary for NoedEdge class objectes (edge) 
    g) nodp : Manager dictionary for NoedEdge class objectes (node) 
    h) flavourp : Manager shared memory object (list) with the free to .pop(0) 
               flavour clolours there is a max flavour set in beginning
    i) fieldp : Manager shared memory object (array) the flavourfield filled 
             elementwise with references type(int) to the flavour of
             NodeEdge class objects. Each edge position is flavour marked in
             the array. The array is a flatten version 
             (1D reshape(dimx,dimy,dimz)) of the field.     
Output:

Inline output:
Plot output:
Save file:
"""
def simple_growth(cod_obj, force_drift, force_struc, view,
                  new_cod, new_nod, new_edd,
                  eddp, nodp, ed_namesp, no_namesp, co_namesp,
                  flavourp, fieldp, frame):
    
    'Draw next growth direction'
    new_posvec = Pdf_Convolver.pdf_convolver(cod_obj, force_drift, force_struc)
    cod_obj.growing(new_posvec)
    
    'update current branching state and time stemp'
    cod_obj.last_splitting +=1
    
    'extent current edge'
    c_edge = cod_obj.current_edge[0]
    new_edd[c_edge]=eddp[c_edge]
    new_edd[c_edge].pos_list=np.concatenate((new_edd[c_edge].pos_list,[cod_obj.pos_list[-1]]),axis=0)


    new_edd[c_edge].x_pos2 = new_edd[c_edge].pos_list[-1][0] 
    new_edd[c_edge].y_pos2 = new_edd[c_edge].pos_list[-1][1] 
    new_edd[c_edge].z_pos2 = new_edd[c_edge].pos_list[-1][2]
    
    'shift cone head(the last node to new growth position)'
    l_node = cod_obj.node_list[-1]
    
    new_nod[l_node] = nodp[l_node]
    new_nod[l_node].pos_list =np.array([cod_obj.pos_list[-1]])
    new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
    new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
    new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2] 
    
    'update eddp proxy'
    eddp.update(new_edd)
    'update nodp proxy'
    nodp.update(new_nod)
    
    'Place edge flavour in flavourfield'
    l_pos = np.round(eddp[c_edge].pos_list[-1])
    dim = cod_obj.field_dim
    x0 = int(l_pos[0])*int(dim[1]*dim[2])
    y0 = int(l_pos[1])*int(dim[2])
    z0 = int(l_pos[2])
    index =  int(x0+y0+z0)     
    fieldp[index] = int(new_edd[c_edge].name[1:])    


"""
Function name : branching_growth()
***Description***

--plain text ---

***I/O***
Input parameter: 
    a) cod_obj : the GrowthCone Class object
    b) force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    c) force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 structure
                  componentes of the substrate
    d) view : a np.array().shape(3,3,3) with information about the next NodeEdge 
           class neighbours
    e) new_cod, new_edd, new_nod : empty dictionaries to intermediate collect 
                                newly constructed GrowthCone class objects and 
                                NodeEdge class objects
    f) eddp : Manager dictionary for NoedEdge class objectes (edge) 
    g) nodp : Manager dictionary for NoedEdge class objectes (node) 
    h) flavourp : Manager shared memory object (list) with the free to .pop(0) 
               flavour clolours there is a max flavour set in beginning
    i) fieldp : Manager shared memory object (array) the flavourfield filled 
             elementwise with references type(int) to the flavour of
             NodeEdge class objects. Each edge position is flavour marked in
             the array. The array is a flatten version 
             (1D reshape(dimx,dimy,dimz)) of the field.
Output:

Inline output:
Plot output:
Save file:
"""
def branching_growth(cod_obj, force_drift, force_struc, view,
                     new_cod, new_nod, new_edd,
                     eddp, nodp, ed_namesp, no_namesp, co_namesp,
                     flavourp, fieldp, frame):
    
    'Call the brachning operation for the cone'
    a, new_cod, new_nod, new_edd = Decision_Maker.branching_decision(cod_obj, force_drift, force_struc, view,
                                                                     new_cod, new_nod, new_edd,
                                                                     eddp, nodp, ed_namesp, no_namesp, co_namesp,
                                                                     flavourp)
    
    'Draw next growth direction for mother cone'
    new_posvec1 = Pdf_Convolver.pdf_convolver(cod_obj, force_drift, force_struc)
    cod_obj.growing(new_posvec1)
    'update current branching state and time stemp'
    cod_obj.last_splitting += 1   
    'extent current edge of the mother cone'

    c_edge = cod_obj.current_edge[0]

    new_edd[c_edge].pos_list=np.concatenate((new_edd[c_edge].pos_list,[cod_obj.pos_list[-1]]),axis=0)


    new_edd[c_edge].x_pos2 = new_edd[c_edge].pos_list[-1][0] 
    new_edd[c_edge].y_pos2 = new_edd[c_edge].pos_list[-1][1] 
    new_edd[c_edge].z_pos2 = new_edd[c_edge].pos_list[-1][2]
    
    'shift cone head(the last node) of the mother cone to new growth position'
    l_node = cod_obj.node_list[-1]
    
    new_nod[l_node].pos_list =np.array([cod_obj.pos_list[-1]])
    new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
    new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
    new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2] 
    

    
    'Place edge of the mother cone flavour in flavourfield'
    l_pos = np.round(new_edd[c_edge].pos_list[-1])
    dim = cod_obj.field_dim
    x0 = int(l_pos[0])*int(dim[1]*dim[2])
    y0 = int(l_pos[1])*int(dim[2])
    z0 = int(l_pos[2])
    index1 =  int(x0+y0+z0)         

    'Draw next growth direction for daughter cone'
    new_posvec2 = Pdf_Convolver.pdf_convolver(new_cod[a], force_drift, force_struc)
    new_cod[a].growing(new_posvec2)
    new_cod[a].frame_seq.append(frame)
    'update current branching state and time stemp'
    new_cod[a].last_splitting += 1
    'extent current edges'
    c_edge2 = new_cod[a].current_edge[0]
    
    new_edd[c_edge2].pos_list=np.concatenate((new_edd[c_edge].pos_list,[new_cod[a].pos_list[-1]]),axis=0)

    
    new_edd[c_edge2].x_pos2 = new_edd[c_edge2].pos_list[-1][0] 
    new_edd[c_edge2].y_pos2 = new_edd[c_edge2].pos_list[-1][1] 
    new_edd[c_edge2].z_pos2 = new_edd[c_edge2].pos_list[-1][2]
    
    'set new node pos for last node // node_list[-1]'
    l_node = new_cod[a].node_list[-1]
    
    new_nod[l_node].pos_list =np.array([new_cod[a].pos_list[-1]])
    new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
    new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
    new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2]          

    'Place edge of the daughter cone flavour in flavourfield'
    l_pos = np.round(new_edd[c_edge2].pos_list[-1])
    dim = cod_obj.field_dim
    x0 = int(l_pos[0])*int(dim[1]*dim[2])
    y0 = int(l_pos[1])*int(dim[2])
    z0 = int(l_pos[2])
    index2 =  int(x0+y0+z0)     
    
    'update eddp proxy'
    eddp.update(new_edd)
    'update nodp proxy'
    nodp.update(new_nod)    

    fieldp[index1] = int(new_edd[c_edge].name[1:])
    fieldp[index2] = int(new_edd[c_edge2].name[1:])  
    


"""
Function name : bifurcation_growth()
***Description***

--plain text ---

***I/O***
Input parameter: 
    a) cod_obj : the GrowthCone Class object
    b) force_drift : a vector, np.array([dx,dy,dz]) with the drift components
    c) force_struc : a vector np.array([xx,xy,xz,yy,yz,zz]) with the 6 structure
                  componentes of the substrate
    d) view : a np.array().shape(3,3,3) with information about the next NodeEdge 
           class neighbours
    e) new_cod, new_edd, new_nod : empty dictionaries to intermediate collect 
                                newly constructed GrowthCone class objects and 
                                NodeEdge class objects
    f) eddp : Manager dictionary for NoedEdge class objectes (edge) 
    g) nodp : Manager dictionary for NoedEdge class objectes (node) 
    h) flavourp : Manager shared memory object (list) with the free to .pop(0) 
               flavour clolours there is a max flavour set in beginning
    i) fieldp : Manager shared memory object (array) the flavourfield filled 
             elementwise with references type(int) to the flavour of
             NodeEdge class objects. Each edge position is flavour marked in
             the array. The array is a flatten version 
             (1D reshape(dimx,dimy,dimz)) of the field.
Output:
	a)
	b)

Inline output:
Plot output:
Save file:
"""
def bifurcation_growth(cod_obj, force_drift, force_struc, view,
                       new_cod, new_nod, new_edd,
                       eddp, nodp, ed_namesp, no_namesp, co_namesp,
                       flavourp, fieldp, frame):
    
    'Call the bifurcation operation for the cone'
    a, b, new_cod, new_nod, new_edd = Decision_Maker.bifurcation_decision(cod_obj, force_drift, force_struc, view,
                                                                         new_cod, new_nod, new_edd,
                                                                         eddp, nodp, ed_namesp, no_namesp, co_namesp,
                                                                         flavourp)
    
    'Draw next growth direction for mother cone'
    new_posvec1 = Pdf_Convolver.pdf_convolver(new_cod[b], force_drift, force_struc)
    new_cod[b].growing(new_posvec1)
    new_cod[b].frame_seq.append(frame)
    'update current branching state and time stemp'
    new_cod[b].last_splitting += 1   
    'extent current edge of the mother cone'

    c_edge1 = new_cod[b].current_edge[0]

    new_edd[c_edge1].pos_list=np.concatenate((new_edd[c_edge1].pos_list,[new_cod[b].pos_list[-1]]),axis=0)


    new_edd[c_edge1].x_pos2 = new_edd[c_edge1].pos_list[-1][0] 
    new_edd[c_edge1].y_pos2 = new_edd[c_edge1].pos_list[-1][1] 
    new_edd[c_edge1].z_pos2 = new_edd[c_edge1].pos_list[-1][2]
    
    'shift cone head(the last node) of the mother cone to new growth position'
    l_node = new_cod[b].node_list[-1]
    
    new_nod[l_node].pos_list =np.array([new_cod[b].pos_list[-1]])
    new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
    new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
    new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2]  

    
    'Place edge of the mother cone flavour in flavourfield'
    l_pos = np.round(new_edd[c_edge1].pos_list[-1])
    dim = new_cod[b].field_dim
    x0 = int(l_pos[0])*int(dim[1]*dim[2])
    y0 = int(l_pos[1])*int(dim[2])
    z0 = int(l_pos[2])
    index1 =  int(x0+y0+z0)     
      

    'Draw next growth direction for daughter cone'
    new_posvec2 = Pdf_Convolver.pdf_convolver(new_cod[a], force_drift, force_struc)
    new_cod[a].growing(new_posvec2)
    new_cod[a].frame_seq.append(frame)
    'update current branching state and time stemp'
    new_cod[a].last_splitting += 1
    'extent current edges'
    c_edge2 = new_cod[a].current_edge[0]
    
    new_edd[c_edge2].pos_list=np.concatenate((new_edd[c_edge2].pos_list,[new_cod[a].pos_list[-1]]),axis=0)

    
    new_edd[c_edge2].x_pos2 = new_edd[c_edge2].pos_list[-1][0] 
    new_edd[c_edge2].y_pos2 = new_edd[c_edge2].pos_list[-1][1] 
    new_edd[c_edge2].z_pos2 = new_edd[c_edge2].pos_list[-1][2]
    
    'set new node pos for last node // node_list[-1]'
    l_node = new_cod[a].node_list[-1]
    
    new_nod[l_node].pos_list =np.array([new_cod[a].pos_list[-1]])
    new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
    new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
    new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2]          

    'Place edge of the daughter cone flavour in flavourfield'
    l_pos = np.round(new_edd[c_edge2].pos_list[-1])
    dim = cod_obj.field_dim
    x0 = int(l_pos[0])*int(dim[1]*dim[2])
    y0 = int(l_pos[1])*int(dim[2])
    z0 = int(l_pos[2])
    index2 =  int(x0+y0+z0)     
    
    'update eddp proxy'
    eddp.update(new_edd)
    'update nodp proxy'
    nodp.update(new_nod)
    
    fieldp[index2]=int(new_edd[c_edge2].name[1:]) 
    fieldp[index1]=int(new_edd[c_edge1].name[1:])  



"""
Function name : merge_edges()
***Description***

This function collects all edges in the veiw of the flavour field and calculates
to one randomly selected edge the shortes distance. A tuple 
as merge_event-variable returns the name of the hit edge, the index of the hit
position on the edge, the current merging cone and the name of the head node
merge_event = (edge_name, index_position, cone_name, node_name)

***I/O***
Input parameter: 
    a) cod_obj : the GrowthCone Class object
    b) view : a np.array().shape(3,3,3) with information about the next NodeEdge 
           class neighbours
    c) eddp : Manager dictionary for NoedEdge class objectes (edge) 
Output:
	a) False if no NodeEdge Class object is near
	b) merge_event : type('tuple') a tuple with all information about the hit
       edge and the hiting position 

Inline output:
Plot output:
Save file:
"""
def merge_edges(cod_obj, view, eddp):
    cone_pos = cod_obj.pos_list[-1]

    nearflavour=view[view > 0]
    nearedge = ['e'+str(j) for j in nearflavour]

    shuffle(nearedge)
    for i in nearedge:
        if eddp[i].flavour != cod_obj.flavour:
            quad_dis = np.sum((eddp[i].pos_list - cone_pos)**2, axis = 1)
                
            if np.amin(quad_dis) <= (cod_obj.growth_length)**2:
                merge_event = [(eddp[i].name, 
                                list(quad_dis).index(np.amin(quad_dis)),
                                cod_obj.name,
                                cod_obj.node_list[-1])]
                return merge_event[0]
    
    return False


"""
Function name : merge_cone_with_edge()
***Description***

If a cone hits another cone or edge its head is moved towards the hit position.

***I/O***
Input parameter: 
    a) cod_obj : the GrowthCone Class object
    b) merge_event : Growthcone Class object attribute type('tuple')
    e) new_cod, new_edd, new_nod : empty dictionaries to intermediate collect 
                                newly constructed GrowthCone class objects and 
                                NodeEdge class objects
    f) eddp : Manager dictionary for NoedEdge class objectes (edge) 
    g) nodp : Manager dictionary for NoedEdge class objectes (node) 
    
Output:

Inline output:
Plot output:
Save file:
"""
def merge_cone_with_edge(cod_obj, merge_event, new_edd, new_nod, eddp, nodp, fieldp):
    'shift  con tip to merge position'
    vector = eddp[merge_event[0]].pos_list[merge_event[1]] - cod_obj.pos_list[-1]
    if np.array_equal(vector, np.array([0,0,0])) == True:
        pass
    else:    
        cod_obj.pos_list =np.concatenate((cod_obj.pos_list,[eddp[merge_event[0]].pos_list[merge_event[1]]]), axis = 0)
        vector = vector / np.linalg.norm(vector)
        cod_obj.vector_list = np.concatenate((cod_obj.vector_list,[vector]), axis = 0)
        
        pol = np.rad2deg(np.arccos(cod_obj.vector_list[-1][2]))
        if pol > 180:
            pol = 360 - pol
        if pol == 0:
            az = 0
        else:
            az = np.rad2deg(np.arctan2(cod_obj.vector_list[-1][1],cod_obj.vector_list[-1][0]))
            if az < 0:
                az = 360 + az
        cod_obj.angle_list = np.concatenate((cod_obj.angle_list,[[pol,az]]),axis = 0)
      
        'extent current edges'
        c_edge = cod_obj.current_edge[0]
        
        new_edd[c_edge]=eddp[c_edge]
        new_edd[c_edge].pos_list=np.concatenate((new_edd[c_edge].pos_list,[cod_obj.pos_list[-1]]),axis=0)
        
        new_edd[c_edge].x_pos2 = new_edd[c_edge].pos_list[-1][0] 
        new_edd[c_edge].y_pos2 = new_edd[c_edge].pos_list[-1][1] 
        new_edd[c_edge].z_pos2 = new_edd[c_edge].pos_list[-1][2]
        
        'set new node pos for last node // node_list[-1]'
        l_node = cod_obj.node_list[-1]
        
        new_nod[l_node] = nodp[l_node]
        new_nod[l_node].pos_list =np.array([cod_obj.pos_list[-1]])
        new_nod[l_node].x_pos = new_nod[l_node].pos_list[-1][0] 
        new_nod[l_node].y_pos = new_nod[l_node].pos_list[-1][1] 
        new_nod[l_node].z_pos = new_nod[l_node].pos_list[-1][2]
        'update eddp proxy'
        eddp.update(new_edd)
        'update nodp proxy'
        nodp.update(new_nod)
        
        'Place edge flavour in flavourfield'
        l_pos = eddp[c_edge].pos_list[-1]
        dim = cod_obj.field_dim
        x0 = int(l_pos[0])*int(dim[1]*dim[2])
        y0 = int(l_pos[1])*int(dim[2])
        z0 = int(l_pos[2])
        index =  int(x0+y0+z0)     
        fieldp[index] = int(new_edd[c_edge].name[1:])
"""
*** THE FOLLOWING FUNCTIONS EXSITS ONLY FOR TESTING AND ARE NOT USEFULL ***
"""
"""
Those following functions are only used in the _main_ of this python file where
growth steps a simulated and tested. The functions are described in the
Initialize_Network_Growth python file. 
"""
def init_manager_init(edd, nod, flavour, field):
    mgr = multiprocessing.Manager()
    eddp = mgr.dict()
    for ke,ve in edd.items():
        eddp[ke]=ve
    nodp = mgr.dict()
    for kn,vn in nod.items():
        nodp[kn]=vn
    flavourp = mgr.list(flavour)
    fieldp = mgr.Array('h',field.reshape(field.shape[0]*field.shape[1]*field.shape[2]))

    return mgr, eddp, nodp, flavourp, fieldp

def external_forces(cod_obj):
    
    'actual cone position'
    x_pos,y_pos,z_pos =cod_obj.pos_list[-1]
    
    'Drift vector components'
    'Drift signal//all directions'
    divx =  divxp([x_pos,y_pos,z_pos])[0]
    divy =  divyp([x_pos,y_pos,z_pos])[0]
    divz =  divzp([x_pos,y_pos,z_pos])[0]

    'Structure tensor components'
    'Diffusion metric // Diffusive tensor components'
    divxx = stpxx([x_pos,y_pos,z_pos])[0]
    divxy = stpxy([x_pos,y_pos,z_pos])[0]
    divxz = stpxz([x_pos,y_pos,z_pos])[0]
    divyy = stpyy([x_pos,y_pos,z_pos])[0]
    divyz = stpyz([x_pos,y_pos,z_pos])[0]
    divzz = stpzz([x_pos,y_pos,z_pos])[0]

    return [divx, divy, divz], [divxx, divxy, divxz, divyy, divyz, divzz]   
    
def view_flavor_field(cod_obj, dim):
    
    l_pos = cod_obj.pos_list[-1]
    "The flavor field is one dimensional array -> the x,y,z pos is flatten"
    x1 = int(l_pos[0]-1)*int((dim[1]*dim[2]))
    x0 = int(l_pos[0])*int((dim[1]*dim[2]))
    x2 = int(l_pos[0]+1)*int((dim[1]*dim[2]))
    y1 = int(l_pos[1]-1)*int((dim[2]))
    y0 = int(l_pos[1])*int((dim[2]))
    y2 = int(l_pos[1]+1)*int((dim[2]))
    z1= int(l_pos[2]-1)
    z0= int(l_pos[2])
    z2= int(l_pos[2]+1)
    
    view =np.array([[[fieldp[x1+y1+z1],fieldp[x0+y1+z1],fieldp[x2+y1+z1]]
                    ,[fieldp[x1+y0+z1],fieldp[x0+y0+z1],fieldp[x2+y0+z1]]
                    ,[fieldp[x1+y2+z1],fieldp[x0+y2+z1],fieldp[x2+y2+z1]]]
                    ,[[fieldp[x1+y1+z0],fieldp[x0+y1+z0],fieldp[x2+y1+z0]]
                    ,[fieldp[x1+y0+z0],fieldp[x0+y0+z0],fieldp[x2+y0+z0]]
                    ,[fieldp[x1+y2+z0],fieldp[x0+y2+z0],fieldp[x2+y2+z0]]]
                    ,[[fieldp[x1+y1+z2],fieldp[x0+y1+z2],fieldp[x2+y1+z2]]
                    ,[fieldp[x1+y0+z2],fieldp[x0+y0+z2],fieldp[x2+y0+z2]]
                    ,[fieldp[x1+y2+z2],fieldp[x0+y2+z2],fieldp[x2+y2+z2]]]])
    return view

if __name__=='__main__':
    
    path_internal = '~/pythrahyper_net/Parameter_Files/internal_parameters.csv'
    path_growth = '~/pythrahyper_net/Parameter_Files/growth_parameters.csv'
    path_multilayer = '~/pythrahyper_net/Parameter_Files/multilayer_dir_parameter.csv'
    path_startpar = '~/pythrahyper_net/Parameter_Files/starting_positions.csv'
    
    import Start_Conditions.Parameter_Importer as Parameter_Importer
    start_pos, st_angle = Parameter_Importer.import_startpos_par(path_startpar)
    internal = Parameter_Importer.import_internal_par(path_internal)
    growth = Parameter_Importer.import_growth_par(path_growth)
    multilayer = Parameter_Importer.import_multilayer_par(path_multilayer)
    
    import Start_Conditions.Node_Maker as Node_Maker 
    flavour= list( np.arange(1,growth['flavour']+1,dtype = int))
    cod = Node_Maker.cone_construction(start_pos, st_angle, internal, growth, multilayer, flavour)
    nod = Node_Maker.node_construction(start_pos)
    cod, nod, edd = Node_Maker.cod_nod_edd_match(cod, nod)
    dim = multilayer['dim']
    field = np.zeros((int(dim[0]),int(dim[1]),int(dim[2])),dtype = int)    
    mgr, eddp, nodp, flavourp, fieldp  =init_manager_init(edd, nod, flavour, field)
    
    import tifffile as tif
    import Start_Conditions.StructField_Maker as StructField_Maker
    imgdx = tif.imread('~/pythrahyper_net/Struct_Field_Images/driftfieldx.tif')
    divxp = StructField_Maker.inter_pol(imgdx,dim)
    imgdy = tif.imread('~/pythrahyper_net/Struct_Field_Images/driftfieldy.tif')
    divyp = StructField_Maker.inter_pol(imgdy,dim)
    imgdz = tif.imread('~/pythrahyper_net/Struct_Field_Images/driftfieldz.tif')
    divzp = StructField_Maker.inter_pol(imgdz,dim)
    imgtxx = tif.imread('~/pythrahyper_net/Struct_Field_Images/strucxx.tif')
    stpxx = StructField_Maker.inter_pol(imgtxx,dim)
    imgtxy = tif.imread('~/pythrahyper_net/Struct_Field_Images/strucxy.tif')
    stpxy = StructField_Maker.inter_pol(imgtxy,dim)
    imgtxz = tif.imread('~/pythrahyper_net/Struct_Field_Images/strucxz.tif')
    stpxz = StructField_Maker.inter_pol(imgtxz,dim)
    imgtyy = tif.imread('~/pythrahyper_net/Struct_Field_Images/strucyy.tif')
    stpyy = StructField_Maker.inter_pol(imgtyy,dim)
    imgtyz = tif.imread('~/pythrahyper_net/Struct_Field_Images/strucyz.tif')
    stpyz = StructField_Maker.inter_pol(imgtyz,dim)
    imgtzz = tif.imread('~/pythrahyper_net/Struct_Field_Images/struczz.tif')
    stpzz = StructField_Maker.inter_pol(imgtzz,dim)
    
    cod_wrap = cod
    for k in range(int(growth['steps'])): 
        new_cod = {}
        new_nod = {}
        new_edd = {}
        print(k)
        "run growth decision, postion recalls and mergings"
        "for each growth cone in dictionary/cluster"
        for key in cod_wrap.keys():
            "Calculate for each growth cone in the ensamble..."
            "..the external forces and the view in the flavor field"
            forces = external_forces(cod_wrap[key])
            view = view_flavor_field(cod_wrap[key], cod_wrap[key].field_dim)
            "Start the growth process"
            cod_obj, new_cod, new_nod, new_edd = growing_function(cod_wrap[key],  #Cone object
                                                                               forces[0],       #Drift
                                                                               forces[1],       #Structure Tensor
                                                                               view,            #View im Flavor Field
                                                                               new_cod,         #Empty Cone dict for actualization
                                                                               new_nod,         #Empty Node dict for actualization
                                                                               new_edd,         #Empty Edge dict for actualization
                                                                               eddp, nodp, flavourp, fieldp, #pointer to the memory dicts/fields
                                                                               k)               #Actual Growth step later used for instance
            cod[key] = cod_obj
        for key in new_cod.keys():
            cod[key] = new_cod[key]
    
        del new_cod
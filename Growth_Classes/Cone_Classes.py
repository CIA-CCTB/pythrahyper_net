# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 17:22:01 2017

@author: Torsten Paul
"""


import numpy as np
from scipy.interpolate import splprep, splev


"""
Class name : ConeClass
***Description***

The main ConeClass class, all further GrowthCone class objects are subclasses
of the class.
IMPORTANT the main class has no __init__ function, it only provides common ConeClass
functions for growing, visualisation, branching and bifurcation.
This class also helds the __del__ function to remove the class object.

Inbuild function:
	- growing the cone move to its next postion
	- __del__ remove class object
	- spline_fitq fits a spline towards the position array
	- spline_visualq creates a 3D plot for the spline
	- branching create the name for the new branching cone
	- bifuration creates the names for both splitting bifuraction
	cones
	- node construction create the name for a new node
	- edge constuction creates the name for a new edge

***I/O***
Input parameter: 
Output:

Inline output:
Plot output:
Save file:
"""
class ConeClass:
    cone_counter = 0
    'destroy as growthcone'
    """
    Function name : __del__()
    ***Description***
    
    Removes the GrowthCone class object the pythonic way for a destructor
    of object orientated programming
    
    ***I/O***
    Input parameter: 
    	a) name the name of the GrowthCone class object
    Output:
    
    Inline output:
    Plot output:
    Save file:
    """  
    def __del__(self):
        ConeClass.cone_counter -= 1

    """
    Function name : growing()
    ***Description***
    
    Class function that performs one growth step for GrowthCone class object.
    The internal pos_list, angle_list and vector_list arrays are extented by the
    imput array selectivly extented by the transformed imput array
    (x,y,z)/norm(x,y,z) -> new_vec
    (x,y,z)/norm(x,y,z) -> (pol,az) -> new_angle
    old_pos +(x,y,z) -> new_pos
    
    ***I/O***
    Input parameter: 
    	a) new_posvec type('np.ndarray').shape(3,) the vector towards the next growth position
	normally a uint vector
    Output:

    Inline output:
    Plot output:
    Save file:
    """  
    def growing(self, new_posvec):

        """
	new_posvec is the new direction in which the cone will grow
	It must be of shape (3,) np.array([x,y,z]) and it is
        transformed in shape (1,3) the growth_vector
	"""
        growth_vector = np.array([new_posvec])
        
        """the polar and azimutal angle are calculated pol [0-pi] az [0-2pi]"""
        pol = np.rad2deg(np.arccos(growth_vector[-1][2]))
        if pol > 180:
            pol = 360 - pol
        if pol == 0:
            az = 0
        else:
            az = np.rad2deg(np.arctan2(growth_vector[-1][1],growth_vector[-1][0]))
            if az < 0:
                az = 360 + az        

        """the polar and azimutal orientation are contected to the angle_list"""
        new_angle = np.array([[pol,az]])
        self.angle_list = np.concatenate((self.angle_list,new_angle),axis = 0)

        """for one growing step an new element is added to the vector_list"""
        self.vector_list = np.concatenate((self.vector_list, growth_vector), axis=0)
        
        
        """superpostion of the last 3 orientations this iskind of an memory for the growing direction"""
        self.vec_mem = np.sum(self.vector_list[-int(self.memory):],axis = 0)/np.linalg.norm(np.sum(self.vector_list[-int(self.memory):],axis = 0))
        
        """
	for one growing step the new position of the the growthcone is added the pos_list
        the orientation is the lastelementrow of the vector_list
	"""
        self.pos_list = np.concatenate(
            (self.pos_list, np.array([self.pos_list[-1]]) + growth_vector * self.growth_length),
            axis=0)
	
    
    """
    Function name : spline_fitq()
    ***Description***
    
    The function calls a fit function splprep from scipy.interpolate. 
    The x,y,z components for each pol_list entry are stored in lists
    xq, yq, zq and transfered to the spline function. A new class attribute
    new_pointsq contains the interplotion parameters to recall the spline
    polynome.
    
    ***I/O***
    Input parameter: 
    Output:

    Inline output:
    Plot output:
    Save file:
    """  
    def spline_fitq(self):

        """" spline fit for the quaternions """
        self.xq = self.pos_list[list(range(self.pos_list.shape[0])), [0] * self.pos_list.shape[0]]

        self.yq = self.pos_list[list(range(self.pos_list.shape[0])), [1] * self.pos_list.shape[0]]
        
        self.zq = self.pos_list[list(range(self.pos_list.shape[0])), [2] * self.pos_list.shape[0]]
        
        tck, u = splprep([self.xq, self.yq, self.zq], s=0)
        self.new_pointsq = splev(u, tck)


    """
    Function name : spline_visualq()
    ***Description***
    
    The function creates a 3D plot Axes3D with the mpl_toolkits.mplot3d
    library. The spline attribute new_pointsq is recalled to plot a rendered
    spline in the plot.
    
    ***I/O***
    Input parameter: 
    Output:
    
    Inline output:
    Plot output:
    Save file:
    """  
    def spline_visualq(self):

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        fig = plt.figure()
        ax3d = fig.add_subplot(111, projection='3d')
        ax3d.plot(self.xq, self.yq, self.zq, 'bx')
        ax3d.plot(self.new_pointsq[0], self.new_pointsq[1], self.new_pointsq[2], 'g-')
        fig.show()
        plt.show()
	
	
    """
    Function name : branching()
    ***Description***
    
    The function creates a branching cone name for a branching event
    the branching_number attribute of the branching cone is increased by 1
    
    ***I/O***
    Input parameter: 
    Output:
    	a) a tpye('str') name of the new branching cone 

    Inline output:
    Plot output:
    Save file:
    """  
    def branching(self, cone1):
        self.branching_number += 1
        a = 'c' + str(cone1)
        return a


    """
    Function name : bifurcation()
    ***Description***
    
    The function creates two bifurction cone names for a bifurcation event
    the bifurcation_number attribute of the cone is increased by 1
    
    ***I/O***
    Input parameter: 
    Output:
    	a) a tpye('str') name of the new bifurcation cone
	b) b tpye('str') name of the new bifurcation cone
    
    Inline output:
    Plot output:
    Save file:
    """  
    def bifurcation(self, cone1, cone2):
        self.bifurcation_number += 1
        a = 'c' + str(cone1)
        b = 'c' + str(cone2)
        return a, b


    """Two functions to construct nodes and edges"""
    """
    Function name : node_constructing()
    ***Description***
    
    The function creates a node name.
    
    ***I/O***
    Input parameter: 
    Output:
    	a) node tpye('str') name of the new node
    
    Inline output:
    Plot output:
    Save file:
    """ 

     
    def node_construction(self, node):
        node = 'n' + str(node)    
    
        return node
    """
    Function name : edge_constructing()
    ***Description***
    
   The function creates a new edge.
    
    ***I/O***
    Input parameter: 
    Output:
    	a) edge type('str') name of the new edge
    
    Inline output:
    Plot output:
    Save file:
    """    
    def edge_construction(self, edge):
        edge = 'e' + str(edge)
        return edge


"""
Class name : GrowthCone
***Description***

First subclass to the ConeClass it inherits all growth, visualisation and
branching functions. 
It adds a _init_ function for GrowthCones class objects and its subclasses 
It defines all GrowthCone class object attributes.

***I/O***
Input parameter: 
Output:

Inline output:
Plot output:
Save file:
"""
class GrowthCone(ConeClass):
    """ __init__ initializes a certain grothcone at pos x,y,z"""
    branching_number = 0
    bifurcation_number = 0 	
    def __init__(self, name, x_pos = 0, y_pos = 0, z_pos = 0):
        ConeClass.cone_counter += 1
        "Internal properties"
        self.name = name
        self.growth_length = 1
        self.aperture = 20
        self.covariance_revert = 40
        self.branching_angle = np.array([90, 0])
        self.bifurcation_angle = np.array([90, 0]) 
        "Cone probabilities for splitting, termination, reactivation..."
        "... and Monte Carlo iterations" 
        self.branchingprob = 0.001
        self.bifurcationprob =0.001
        self.terminationprob = 0.0001
        self.reactivationprob = 0.0001
        self.montecarlo_iterations = 100000
        "Starting point of the growthcone"
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos
        self.pos_list = np.array([[self.x_pos, self.y_pos, self.z_pos]])
        self.vector_list = np.array([[0, 1, 0]])
        self.vec_mem = np.array([[0, 1, 0]])
        self.angle_list = np.array([[90,90]])
        "Memory and imformation about nodes, edges and splitting events"
        "and infomation about the compfield/substrate and search field"      
        self.field_dim = np.array([100,100,100])
        self.max_drift = 1
        self.min_drift = 0
        self.max_eig = 1
        self.min_eig = 0        
        self.searchfield = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])               
        self.flavour = 0
        self.node_list = []
        self.current_edge = []
        self.sleep_merge_event= [('growing')]
        self.last_splitting = 0
        self.frame_seq = []
        
        
"""
Class name : BranchingCone()
***Description***

First subclass to GrowthCone class it inherits all growth, visualisation,
branching, bifucation and construction functions. 
It overwrites the _init_ function of the GrowthCones class
It defines all BranchingCone class object attributes.

***I/O***
Input parameter: 
Output:

Inline output:
Plot output:
Save file:
"""       
class BranchingCone(GrowthCone):
    """ __init__ initializes a certain grothcone at pos x,y,z"""

    def __init__(self, name, x_pos = 0, y_pos = 0, z_pos = 0):
        ConeClass.cone_counter += 1
        "Internal properties"
        self.name = name
        self.growth_length = 1
        self.aperture = 20
        self.covariance_revert = 40
        self.branching_angle = np.array([90, 0])
        self.bifurcation_angle = np.array([90, 0]) 
        "Cone probabilities for splitting, termination, reactivation..."
        "... and Monte Carlo iterations" 
        self.branchingprob = 0.001
        self.bifurcationprob =0.001
        self.terminationprob = 0.0001
        self.reactivationprob = 0.0001
        self.montecarlo_iterations = 100000        
        "Starting point of the growthcone"
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos
        self.pos_list = np.array([[self.x_pos, self.y_pos, self.z_pos]])
        self.vector_list = np.array([[0, 1, 0]])
        self.vec_mem = np.array([[0, 1, 0]])
        self.angle_list = np.array([[90,90]])
        "Memory and imformation about nodes, edges and splitting events"
        "and infomation about the compfield/substrate and search field"      
        self.field_dim = np.array([100,100,100])
        self.max_drift = 1
        self.min_drift = 0
        self.max_eig = 1
        self.min_eig = 0
        self.searchfield = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])          
        self.flavour = 0
        self.node_list = []
        self.current_edge = []
        self.sleep_merge_event= [('growing')]
        self.last_splitting = 0
        self.frame_seq = []
  

"""
Class name : BifurcationCone()
***Description***

Second subclass to GrowthCone class it inherits all growth, visualisation,
branching, bifucation and construction functions. 
It overwrites the _init_ function of the GrowthCones class
It defines all BifurationCone class object attributes.

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
class BifurcationCone(GrowthCone):
    """ __init__ initializes a certain grothcone at pos x,y,z"""

    def __init__(self, name, x_pos = 0, y_pos = 0, z_pos = 0):
        ConeClass.cone_counter += 1
        "Internal properties"
        self.name = name
        self.growth_length = 1
        self.aperture = 20
        self.covariance_revert = 40
        self.branching_angle = np.array([90, 0])
        self.bifurcation_angle = np.array([90, 0]) 
        "Cone probabilities for splitting, termination, reactivation..."
        "... and Monte Carlo iterations"
        self.branchingprob = 0.001
        self.bifurcationprob =0.001
        self.terminationprob = 0.0001
        self.reactivationprob = 0.0001
        self.montecarlo_iterations = 100000
        "Starting point of the growthcone"
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos
        self.pos_list = np.array([[self.x_pos, self.y_pos, self.z_pos]])
        self.vector_list = np.array([[0, 1, 0]])
        self.vec_mem = np.array([[0, 1, 0]])
        self.angle_list = np.array([[90,90]])
        "Memory and imformation about nodes, edges and splitting events"
        "and infomation about the compfield/substrate and search field"      
        self.field_dim = np.array([100,100,100])
        self.max_drift = 1
        self.min_drift = 0
        self.max_eig = 1
        self.min_eig = 0
        self.searchfield = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])               
        self.flavour = 0
        self.node_list = []
        self.current_edge = []
        self.sleep_merge_event= [('growing')]
        self.splitting = 0
        self.frame_seq = []

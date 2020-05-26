#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 08:39:19 2019
@author: top40ub
"""
#import numpy
import numpy as np 
import csv
import os
from pathlib import Path
"""
from the Math_and_Simulation folder to .py-files are loaded containing
functions to create a point clouds arround positions in 3D
"""
from Math_and_Simulation import Pdf_Convolver as PdfC
from Math_and_Simulation import Elliptically_symmetric_Angular_G as EsAG

""" 
    Parameter set toconstruct the multilayered field:
The path to the image used as computation field
The path to the image from which the gradient fields are calculated
The path to the image to derive the structure tensor values
"""

image_dir = {'image' : '/Structured_Images/sphereshell2.tif',
            'imaged' : '/Structured_Images/sphereshell2.tif',
            'imaget' : '/Structured_Images/sphereshell2.tif',
            'image_name' : 'volume.tif',
            'dirft_image_name' : 'volume.tif',
            'strucure_image_name' : 'ramp.tif',
            }
""" 
    The set of parameters for the growth cone
The number of growth steps during one instance 
The number of growth instances
The number of possible flavours
The number of possible edges
The number of possible nodes
The number of possible cones
The number of max. Monte Carlo iterations for the creation of new growth 
orientations
The radius of the Delauny Cluster
"""

"""Those are the default key-value pairs of the growth variable dictionary"""
growth_var = {'steps' : 10,
              'instances' :25,
              'flavour' : 100000,
    	      'edge_number': 100000,
    	      'cone_number': 100000,
    	      'node_number': 100000,
              'covariance_revert' : 40,
              'montecarlo_iterations' : 100000,
              'delauny_cluster_radius' : 10,
              'Proxy_drift': False,
              'Proxy_tensor': False,
              'Proxy_corr': False,
              'Proxy_reverse_eig': True
              }

"""
    The set of parameters for the growth cone:
The aperture of the solid angle: persitence characteristic
The step size: hard coded as integers default should be 1
The splitting angles for branching and bifurcation
The searchfield for a cone to check its surroundings for other
network edges
The probabilities to branch, bifurcate, terminate, or reactivate 
The correlation memory, the number of last growth directions to
correlate with the next one
"""

"""Those are the default key-value pairs for the internal growth parameter"""
internal_var = {'aperture' : 30,
                'stepsize' : 1,
                'branching_angle' : np.array([75,0]),
                'bifurcation_angle' : np.array([75,0]),
                'searchfield' : np.array([0., 0., 0.,
                                           0., 0., 10.,
                                           0., 0., 0.,
                                           0., 0., 10.,
                                           0., 10., 10.,
                                           0., 0., 10.,
                                           0., 0., 0.,
                                           0., 0., 10.,
                                           0., 0., 0.]),
                'branchingprob' : 0.02,
                'bifurcationprob' : 0.005,
                'deathprob' : 0.0001,
                'reactivationprob' : 0.0001,
                'memory': 5
                }  


"""
Function name : start_pos_orientation()
***Description***
This function creates the initial cone start positions and their orientations.
An array with shape (n,3) contains the postions of each cone (x,y,z), an
array (n,2) its orientation (pol,az), a last one (n,1) indicates a common 
start node. 
The function calls a startpoint generator that creates a equal distributed point
cloud around the start position. This cloud is highly tuneable is terms of:
	- shape via aperture and offset
	- density per surface element with point number
	and metric of space via the axes scaling
A switch called common_node indicates if all cones of one cloud belong to one
common node
For more information about the called function start_point_generator see its
describtion in Math_and_Simultation.Elliptically_symmetric_Angular_G
***I/O***
Input parameter: 
Output:
	a) start_pos type('ndarray') shape (n,3) with the cone positions
	b) start_angle type('ndarray') shape (n,2) with orientation
	c) common_node type('ndarray') shape (n,1) with the name of the start
	node
Inline output:
Plot output:
Save file:
"""
def start_pos_orientation():
    
    start_angle_list = np.array([np.random.rand(2)*np.array([180,360]) \
                            for i in range(4)])
    start_pos_list = np.array([np.random.rand(3)*np.array([100,100,100]\
                          ) + np.array([78,78,78]) for i in range(4)])
    
    
    start_pos_list = np.array([[125,125,27],
                               [125,125,223],
                               [125,27,125],
                               [125,223,125],
                               [27,125,125],
                               [223,125,125]])
   
    start_angle_list = np.array([[0,0],
                                 [0,0],
                                 [90,90],
                                 [90,90],
                                 [90,0],
                                 [90,0]])
    
    
    start_pos_list = np.array([#[42,42,85],[127,42,85],[212,42,85],
                               [42,127,85],[127,127,85],[212,127,85],
                               #[42,212,85],[127,212,85],[212,212,85],
                               #[42,42,170],[127,42,170],[212,42,170],
                               [42,127,170],[127,127,170],[212,127,170],
                               #[42,212,170],[127,212,170],[212,212,170]
                               ])
    start_angle_list = np.array([#[0,0],[0,0],[0,0],
                                 #[0,0],[0,0],[0,0],
                                 #[0,0],[0,0],[0,0],
                                 #[0,0],[0,0],[0,0],
                                 [0,0],[0,0],[0,0],
                                 [0,0],[0,0],[0,0]])
    
    start_pos_list = np.array([[124,124,124]])
    start_angle_list = np.array([[0,0]])

    start_pos, start_angle,\
    common_node = EsAG.start_point_generator(start_pos_list,
                                            start_angle_list,
                                            aperture = 180, offset = 0,
                                            points_on_sphere = 30, 
                                            metric = np.array([[100,100,100]]),
                                            common_starting_node = False,
                                            tangent_orientation = 'pol')
    return start_pos, start_angle, common_node


"""
Function name : create_starting_csv()
***Description***
Creates a .csv field with name and direction stored in path_paramter.
The function takes the growth_var dictionary as imput to create a .csv
file with the dictionary key, value pairs as rows.
A header describts the content of the .csv file
***I/O***
Input parameter: 
	a) growth_var type('dict') dictionary with the network growth parameter
	b) path_paramter path and filename of the .csv-file to create
Output:
Inline output:
Plot output:
Save file: Save the .csv-file to directory and with name in path_paramter
"""
def create_starting_csv(growth_var, path_parameter):
    with open(path_parameter, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['The growth parameter for the network growth'])
        for key,value in growth_var.items():        
            if type(value) is np.ndarray:
                filewriter.writerow([key, *value])
            else: 
                filewriter.writerow([key, value])


"""
Function name : create_structured_images_dir_csv()
***Description***
Creates a .csv field with name and direction stored in path_image.
The function takes the image_dir dictionary as imput to create a .csv
file with the dictionary key, value pairs as rows.
A header describts the content of the .csv file
***I/O***
Input parameter: 
	a) image_dir type('dict') dictionary with the the path to the
	multilayer images and their names
	b) path_image path and filename of the .csv-file to create
Output:
Inline output:
Plot output:
Save file: Save the .csv-file to directory and with name in path_image
"""
def create_structured_images_dir_csv(image_dir, path_image):
    with open(path_image, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['The directories of the images to construct the multilayered conmputation field'])
        for key_value, i in zip(image_dir.items(), range(len(image_dir))): 
            if type(key_value[1]) is np.ndarray:
                filewriter.writerow([key_value[0], *key_value[1]])
            else:
                if i <= 2:
                    filewriter.writerow([key_value[0], home + key_value[1]])
                else:
                    filewriter.writerow([key_value[0], key_value[1]])


"""
Function name : create_internal_csv()
***Description***
Creates a .csv field with name and direction stored in path_internal.
The function takes the internal_var dictionary as imput to create a .csv
file with the dictionary key, value pairs as rows.
A header describts the content of the .csv file
***I/O***
Input parameter: 
	a) internal_var type('dict') dictionary with the network growth internal
	parameter
	b) path_internal path and filename of the .csv-file to create
Output:
Inline output:
Plot output:
Save file: Save the .csv-file to directory and with name in path_internal
"""
def create_internal_csv(internal_var, path_internal):
    with open(path_internal, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['The internal parameter for the network growth'])
        for key,value in internal_var.items():
            if type(value) is np.ndarray:
                filewriter.writerow([key, *value])
            else: 
                filewriter.writerow([key, value])


"""
Function name : create_starting_pos()
***Description***
Creates a .csv field with name and direction stored in path_pos.
The function takes the internal_var dictionary as imput to create a .csv
file with n rows for n growth cones and 6 columns for (x,y,z,pol,az,node)
the start position, orientation and start node.
A header describts the content of the .csv file
***I/O***
Input parameter: 
	a) start_pos type('ndarray') shape (n,3) with the cone positions
	b) start_angle type('ndarray') shape (n,2) with orientation
	c) common_node type('ndarray') shape (n,1) with the name of the start
	node
	d) path_pos path and filename of the .csv-file to create
Output:
Inline output:
Plot output:
Save file: Save the .csv-file to directory and with name in path_pos
"""          
def create_starting_pos(path_pos, start_pos, start_angle, common_node):
    with open(path_pos, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(["The cones' start positions and orientations for the network growth"])
        filewriter.writerow(['X','Y','Z','Az','Pol', 'CN'])
        
        for i in range(start_pos.shape[0]):
            filewriter.writerow([*start_pos[i],*start_angle[i],*common_node[i]])


if __name__=='__main__':
    home = str(Path(os.getcwd()).parent)

    
    
    """Creates start positions, start directions growthaperture, stepsize andoverall step number"""
    path_image = home + '/Parameter_Files/structured_image_dir.csv'
    path_parameter = home +  '/Parameter_Files/growth_parameters.csv'
    path_internal = home + '/Parameter_Files/internal_parameters.csv'
    path_pos = home + '/Parameter_Files/starting_positions.csv'
    create_starting_csv(growth_var, path_parameter)
    create_starting_pos(path_pos, *start_pos_orientation())
    create_internal_csv(internal_var, path_internal)
    create_structured_images_dir_csv(image_dir, path_image)
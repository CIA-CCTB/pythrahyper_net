#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os

simdir = sys.argv[1]
run = sys.argv[2]

home = os.getcwd();

#basedir = "/Examples/Osteocytes/";
#simdir = "memory/20/";
#rundir = "1/";

import numpy as np
from time import time as timet
import concurrent.futures
import os
import itertools
import multiprocessing
from pathlib import Path

import tifffile as tif


# Import the framework's main module
import Initialize_Network_Growing as ing

# This python file loads the parameter setting for the growth simulation
import Start_Conditions.Parameter_Importer as Parameter_Importer

# This python file initializes the starting nodes and growth cones: 
import Start_Conditions.Node_Maker as Node_Maker

# This loads the images containing the structure of teh environment:
import Start_Conditions.StructField_Maker as StructField_Maker

# Growth step functionality
import Math_and_Simulation.Step_Maker as Step_Maker

# Functions for plotting and making movies
import Plots_and_Movies.Plot_Maker as Plot_Maker
import Plots_and_Movies.Movie_Maker as Movie_Maker

# Functions for generating the networkx graph
import Network_Graph_Construction.Network_X as Network_X
import Network_Graph_Construction.Update_Network as Update_Network


# Path to this notebook
home = os.getcwd() 

# Locations of csv configuration files
path_multilayer = home + '/Parameter_Files/multilayer_dir_parameters.csv'
path_startpar = home + simdir + 'starting_positions_ocy.csv'
path_growth = home + simdir + 'growth_parameters.csv'
path_internal = home + simdir + 'internal_parameters.csv'
path_structured_images = home + simdir + 'structured_image_dir.csv'



# Generate feature maps from image data of growth environment
#features = ing.StructField_Maker.structured_field(path_structured_images, home, sigma_divd=2, sigma_divt1=2, sigma_divt2=2)

# Initialise the computation grid
field_ref = ing.StructField_Maker.interpol_external_fields(path_multilayer)

# Initialize the object dictionaries
obj_par_env = ing.init_objects_parameter_enviroment(path_startpar,path_internal,path_growth,path_multilayer)



# Extract individual dictionaries
cod, nod, edd, flavour, ed_names, no_names, co_names, field, steps, instances, dim, radius  = obj_par_env

# Create shared memory proxy objects for all parameter and class objects
mgr, eddp, nodp, flavourp, ed_namesp, \
no_namesp, co_namesp, fieldp, forcep = ing.init_manager_init(edd, nod, flavour,
                                                            ed_names, no_names,
                                                            co_names, field,
                                                            field_ref)

# Transfer proxy objects to the correct namespace (needed when running in notebook)
ing.mgr = mgr
ing.eddp = eddp
ing.nodp = nodp
ing.flavourp = flavourp
ing.ed_namesp = ed_namesp
ing.no_namesp = no_namesp
ing.co_namesp = co_namesp
ing.fieldp = fieldp
ing.forcep = forcep



# Starting the growth process simulation

#%matplotlib inline

instances = 10
steps = 5

growing_results = ing.init_par_growing(cod, nod, edd, steps,instances, dim, radius,eddp, nodp, flavourp, fieldp, ed_namesp, no_namesp, co_namesp,forcep,
                                  timer = True,
                                  plot = False,
                                  Voronoi = False,
                                  network_x = False)

cod, nod, edd, flavour, ed_names, no_names, co_names, vor, G_vor, subG_vor, G, subG = growing_results




import pickle


rundir = str(run)+"/"

run_path = home + simdir + rundir;

if not os.path.exists(run_path):
    os.makedirs(run_path);
    print(run_path+"generated!")

with open(run_path + 'cod.pkl', 'wb') as f:
    pickle.dump(cod, f)
    
with open(run_path + 'edd.pkl', 'wb') as f:
    pickle.dump(edd, f)
    
with open(run_path + 'nod.pkl', 'wb') as f:
    pickle.dump(nod, f)
    
import networkx as nx

G = nx.Graph()
    
nodes = [key for key in nod.keys()]    

G.add_nodes_from(nodes)
    
edges = [(edd[key].nodes[0],edd[key].nodes[1],len(edd[key].pos_list)) for key in edd.keys()]

G.add_weighted_edges_from(edges)

thr = 1;
for key in edd.keys():
    dx = edd[key].x_pos1 - edd[key].x_pos2;
    dy = edd[key].y_pos1 - edd[key].y_pos2;
    dz = edd[key].z_pos1 - edd[key].z_pos2;
    d = np.sqrt(dx**2 + dy**2 + dz**2);
    n1,n2 = edd[key].nodes[:]        
    if (d<thr) and (n1 in G) and (n2 in G):
        G = nx.contracted_nodes(G,edd[key].nodes[0],edd[key].nodes[1],);
    
G.remove_edges_from(nx.selfloop_edges(G))
        
with open(run_path + 'G.pkl', 'wb') as f:
    nx.write_gpickle(G, f)

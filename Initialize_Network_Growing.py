# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 15:16:56 2017

@author: Torsten Paul
"""

"""THIS IS THE INITIZALISATION FILE FOR THE WHOLE GROWTH PROCESS"""

import numpy as np
from time import time as timet
import concurrent.futures
import os
import itertools
import multiprocessing
from pathlib import Path

"""The python files containing all needed self writen function"""

"""This python file initializes the starting nodes and with their starting 
direction an creates the first cone dictonary"""
import Start_Conditions.Node_Maker as Node_Maker
"""This python file let the cone grow, also there are occur decisions for 
branching and bifucation together with feed back from the enviroment"""
import Math_and_Simulation.Step_Maker as Step_Maker
"""This python file loads the image containing the enviroment structure and 
creates the computational field with the struc tensor representation"""
import Start_Conditions.StructField_Maker as StructField_Maker
"""After the growth process the cones and their path ways are plotted"""
import Plots_and_Movies.Plot_Maker as Plot_Maker
import Plots_and_Movies.Movie_Maker as Movie_Maker
"""This python file loads the parameter setting for the growth simulation"""
import Start_Conditions.Parameter_Importer as Parameter_Importer
"""These python files contain functions to construct the mathematical 
network graph"""
import Network_Graph_Construction.Network_X as Network_X
import Network_Graph_Construction.Update_Network as Update_Network


"""
Function name : init_objects_parameter_enviroment()
***Description***

Transfers the parameters from all parameter files (<- input path strings) in 
dictionaries or lists (-> start_pos, st_ angle, internal, growth, multilayer)
From those parameter dictionaries the class object dictionaries cod, nod, edd
and the growth simulation objects/parameter flavour, dim, steps, instances,
computation field and the cluster radius are initialised.

flavour (<- growth)
cod (<- start_pos, start_pos, st_ angle, internal, growth, multilayer, flavour)
nod (<- start_pos)
edd (<- cod, nod)
Each GrowthCone class object has its own flavour colour and places a reference
to it at its starting position in the computation field, the flavourfield

***I/O***
Input parameter: 
	a) path_startpar type('str') path to the starting_position.csv file
	b) path_internal type('str') path to the internal_paramters.csv file
	c) path_growth type('str') path to the growth_paramters.csv file
	d) path_multilayer tpe('str') path to the multilayer_dir_csv file
Output:
	a) cod type('dict') dictionary that contains the initial GrowthCone 
	class objects keys:'c0',..,'cn'
	b) nod type('dct') dictionary that contains the initial NodeEdge class
	objects (node) 'n0','n0-1',.., 'nn','nn-1' 
	c) edd type('dict') dictionary that contains the initial NodeEdge class
	objects (edge) keys: '1',..,'n+1'
	d) flavour type('list') a list contains the range of flavour colours
	e) edge_names type('list') a list contains the range of edge names
	f) node_names type('list') a list contains the range of node_names
	g) cone_names type('list') a list contains the range of cone_names
	h) field type('ndarray().shape(dim)') the computation field
	i) steps type('array()') the number of growth iterations per instance
	j) instances type('array()') the number of instance repititions
	k) dim type('ndarray().shape(3)) the dimension of the computation field
	l) radius_clus type('array()') the radius for the delauny graph cluster

Inline output: print message ('number of starting cones')
Plot output:
Save file:
"""
def init_objects_parameter_enviroment(path_startpar, 
                                      path_internal, 
                                      path_growth, 
                                      path_multilayer):

    start_pos, st_angle, common_node = Parameter_Importer.import_startpos_par(path_startpar)
    internal = Parameter_Importer.import_internal_par(path_internal)
    growth = Parameter_Importer.import_growth_par(path_growth)
    multilayer = Parameter_Importer.import_multilayer_par(path_multilayer)
    
    flavour= list(np.arange(1,growth['flavour']+1,dtype = int))
    edge_names = list(np.arange(1,growth['edge_number']+1,dtype = int))
    node_names = list(np.arange(1,growth['node_number']+1,dtype = int))
    cone_names = list(np.arange(1,growth['cone_number']+1,dtype = int))
    cod = Node_Maker.cone_construction(start_pos, st_angle, 
                                       internal, growth, multilayer,
                                       flavour, cone_names)
				       
    
    nod = Node_Maker.node_construction(start_pos, common_node, node_names)
    cod, nod, edd = Node_Maker.cod_nod_edd_match(cod, nod, common_node, 
                                                 node_names, edge_names)
    dim = multilayer['dim']
    steps = int(growth['steps'][0])
    instances = int(growth['instances'][0])
    radius_clus = growth['delauny_cluster_radius']
    field = np.zeros((int(dim[0]),int(dim[1]),int(dim[2])),dtype = int)
    
    """Place cone flavour in flavourfield"""
    for key in edd:
        cone_pos = edd[key].pos_list[-1]
        x0 = int(cone_pos[0])
        y0 = int(cone_pos[1])
        z0 = int(cone_pos[2])
        field[x0,y0,z0] = int(edd[key].name[1:])
    print('The number of starting cones is :',len(cod))
    return cod, nod, edd, flavour, edge_names, node_names, \
           cone_names, field, steps, instances, dim, radius_clus


"""
Function name : dict_wrapper_growth()
***Description***

The dict_wrapper_growth() function is an envelop for the growing_function().
It takes as input a key string and the corresponding GrowthCone class object as
well as the step iteration parameter. It reconstructs a new dictionary cod_wrap
for the key and GrowthCone object. 
1: Loop runs over the number of iteration steps.
    In each iteration three empts auxiliary dictionaries new_cod, new_nod, 
    new_edd are initialised.
    2: A nested loop runs over each key/GrowthCone pair in cod_wrap.
    First checks if the GrowthCone class objcet is still able to grow
    In side the second loop each growth cone calls its surrounding forces and
    takes a view of the surrounding network structure.
    The function growing_function() (<- cod_wrap, forces, view, new_cod,
                                     new_nod, new_edd, shared memory objects)
    Retruns the GrowthCone class object in a new state/position
    if new GrowthCone class objects or NodeEdge Class objects are created they
    are stored in the auxiliary dictionaries. New GrowthCone class objects are
    transfered to cod_wrap
Returns the cod_warp dictionary after the loops    

***I/O***
Input parameter: 
	a) *args: args[0] a cod key type('str')
               args[1] a cod[key] value a GrowthCone class objects
               args[2] steps type('array()') the number of growth iterations
               args[3] integer with value of the current instance
       further args[4:] are optional
	b) **kwargs empty
Output:
	a) cod_wrap type('dict') dictionary that contains {arg[0]:arg[1]} and 
        further constructed GrowthCone class objects

Inline output:
Plot output:
Save file:
"""
def dict_wrapper_growth(*args,**kwargs):
    cod_wrap = {args[0]:args[1]}
    instance = int(args[3])
    steps = int(args[2])
    for k in range(steps):
        """new auxiliary dictionaries for branching/bifuration cones"""
        """new nodes and new edges"""
        new_cod = {}
        new_nod = {}
        new_edd = {}
        """run growth decision, postion recalls and mergings"""
        """for each growth cone in dictionary/cluster"""
        for key in cod_wrap.keys():
            """Check if cone is able to grow"""
            if cod_wrap[key].sleep_merge_event == [('growing')] or \
                cod_wrap[key].sleep_merge_event == [('sleeping')]:
                """Calculate for each growth cone in the ensamble..."""
                """..the external forces and the view in the flavor field"""
                forces = external_forces(cod_wrap[key])
                view = view_flavour_field(cod_wrap[key],
                                          cod_wrap[key].field_dim)
                """Start the growth process"""
                frame = k + instance * steps 
                class_obj = Step_Maker.growing_function(cod_wrap[key],  
                                                        forces[0], forces[1],
                                                        view, new_cod, new_nod,
                                                        new_edd, eddp, nodp,
                                                        ed_namesp, no_namesp, co_namesp,
                                                        flavourp, fieldp, frame)
                cod_obj, new_cod, new_nod, new_edd = class_obj
        for key in new_cod.keys():
            cod_wrap[key] = new_cod[key]
    return cod_wrap


"""
Function name : build_arg_map()
***Description***

An auxiliary function to help entangle the cod object as largest dictionary and
all other parameter from *arg arguments.
The function build_arg_map builds a map/list (arg_map) containing zeros for 
each element in args. 
Loop over all elements in args[:]:
    If the ith element is type('dict')
    Assign the length len(args[i]) to the map/list (arg_map) object's
    zeroth element [len(dict),0,...,0]

Args should only have one element with type('dict'), if otherwise the cod
object should be the largest dictionary.

max(arg_map) is later used as the number processes and hard copys of each 
args[1:] element to pass to each specific process. 

***I/O***
Input parameter: 
	a) *args: args[0] cod tpye('dict) containing the GrowthCone class objects
               args[1] steps type('array()') the number of growth iterations
       further args[2:] are optional
	b) **kwargs empty
Output:
	a) arg_map type('list') [len(dict),0,...,0]

Inline output:
Plot output:
Save file:
"""
def build_arg_map(*args):
    arg_map = [0 for i in range(len(args))]
    for e in enumerate(args):
        
        if type(e[1]) == dict:
            arg_map[e[0]] = len(e[1])
        else:
            pass
        
    return arg_map


"""
Function name : multiprocessor_dictp()
***Description***

The multiprocessor_dictp() takes the parameter cod, steps and other optional 
parameter as arguments. If arguments other than the cod dictionary are of 
type('dict') len(cod) needs to be the largest dictionary.
The functions initialises a multiprocessing excecution of the function 
dict_wrapper_growth. It sets all avaliable CPUs as process workers.

The function build_arg_map runs over all args and constructs 2 iter_arg lists
iter_arg1 and iter_arg2. Both have the same lenght len(max(arg_map)) and 
contain the keys and values(GrowthCone class objects) of cod. 

The concurrent PoolExcecutorh key/value pair to the worker pool.
The function dict_wrapper_growth creates a single dictionary with only one 
key/value pair from each pooled key/value pair of cod together with a hard copy
of each other argument the growing_function is excetuted.

***I/O***
Input parameter: 
	a) *args: args[0] cod tpye('dict) containing the GrowthCone class objects
               args[1] steps type('array()') the number of growth iterations
       further args[2:] are optional
	b) **kwargs empty
Output:
	a) list(exlist) type('list') list with dictionaries of GrowthCone class 
        objects as elements.

Inline output:
Plot output:
Save file:
"""
def multiprocessor_dictp(*args, **kwargs):
    WORKERS = os.cpu_count()
    WORKERS = 4
    arg_map = build_arg_map(*args)
    repeats = max(arg_map)
    

   
    iter_arg1 = list(args[arg_map.index(max(arg_map))].keys())
    iter_arg2 = list(args[arg_map.index(max(arg_map))].values())
    
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=WORKERS) as exe:
        exlist = exe.map(dict_wrapper_growth, 
                         iter_arg1, 
                         iter_arg2, 
                         *[itertools.repeat(element, repeats) \
                           for element in args[1:]])

    return list(exlist)


"""
Function name : init_par_growing()
***Description***

The function starts the multiprocessing network growth. Running a
loop up to instances the function multiprocessor_dictp is excecuted with
the cone dictionary cod and the number of steps as variables.
Eddp and nodp as well flavourp and fieldp are proxy manager objects to whom
each process has access. These shared memories are accompanied by the features
of the compfield:
	the gradients div**
	the structure tensor stp** 
Those are function operators interpolation over the multilayer feature field
stored in the common namespace.

After all growth instances the NodeEdge class objects are updated and matched
with the proxy dictionaries. There are four switches to further elucidate the
network structure.

Four switches tpye('Boolian'):
    plot: plot the network morphology in a 3D plot
    Voronoi: Clusters all active growth cones within the possible contact
    radius of value radius_clus  
    timer: to measure the run time and prints it
    network_x: activates the interface to the Network_X library and transfers
    the NodeEdge Class objects into a Network_X graph.

***I/O***
Input parameter: 
	a) cod type('dict') dictionary that contains the initial GrowthCone 
	class objects keys:'c0',..,'cn'
	b) nod type('dct') dictionary that contains the initial NodeEdge class
	objects (node) 'n0','n0-1',.., 'nn','nn-1' 
	c) edd type('dict') dictionary that contains the initial NodeEdge class
	objects (edge) keys: '1',..,'n'
	d) steps type('array()') the number of growth iterations per instance
	e) instances type('array()') the number of instance repititions
	f) dim type('ndarray().shape(3)) the dimension of the computation field
	g) radius_clus type('array()') the radius for the delauny graph cluster
	h) eddp proxy type('dict') of edd
	i) nodp proxy type('dict') of nod
	j) flavourp proxy type('list') of flavour list
	k) fieldp proxy type('array) 1D (dim[0]*dim[1]*dim[2]) of field
	l) ed_namesp proxy type('list') a list contains the range of edge names
	m) no_namesp proxy type('list') a list contains the range of node_names
	n) co_namesp proxy type('list') a list contains the range of cone_names	
	o) stpxx, stpxy, stpxz, stpyy, stpyz, stpzz function operators for the parts
	of the multilayer structure tensor
	p) divxp, divyp, divzp function operators for the parts of the multilayer 
	drift field
	q) four switches type('Boolean') activate or inactivate timing, plotting, 
	Voronoi clustering and network analysis
	
Output:
	a) cod type('dict') dictionary that contains all exsisting GrowthCone 
	class objects
	b) nod type('dct') dictionary that contains all exsisting NodeEdge class
	objects (node)
	c) edd type('dict') dictionary that contains all exsiting NodeEdge class
	objects (edge)
	d) flavour type('list') reduced list where flavours hs been distributed via
	.pop(0) to new NodeEdge class objects (edge)
	e) ed_names type('list') reduced list contains the range of edge names
	f) no_names type('list') reduced list contains the range of node_names
	g) co_names type('list') reducec list contains the range of cone_names
    h) vor, G_vor, subG_vor type('list') empty if Voronoi switch is False, 
    Voronoi cluster if True
    i) G, subG type('list') empty if Networ_X switch is False, network_X graph
    if True

Inline output: Prints the number for starting growth cones per instances and number 
of growth steps per instance, running time if swtchis active, after graph
construction the resulting edge and node number
Plot output: Plots the 3D network, the Voronoi graph per instance and the network
topology graph if the particular switches are active
Save file:
"""
def init_par_growing(cod, nod, edd, steps, instances, dim, radius,
                     eddp, nodp, flavourp, fieldp,
                     ed_namesp, no_namesp, co_namesp,
                     forcep, timer = False, plot = False, 
                     Voronoi = False, network_x = True ):
    
    tic = timet()
    for i in range(int(instances)):
        print('Instance number ', i, ' with ', steps, ' steps', ' and ',
              len(cod), 'cones, ',
              len(eddp), 'edges and ',
              len(nodp), 'nodes.')
        """Executing of the parallel network growth"""

        """
        calls the step_maker file 
        and starts a multiprocessing growth process
        """
        podp = multiprocessor_dictp(cod, steps, i)
        """after all parallel steps the cod is updated"""
        for i in podp:
            cod.update(i)

    tac = timet()
    
    """eddp, and nodp are manager dicts which could be read and written by each 
    pool worker. There key/value pairs are updated to edd and nod"""
    for key in eddp.keys():
        edd[key]=eddp[key]
    for key in nodp.keys():
        nod[key]=nodp[key]
    flavour = [i for i in flavourp]
    ed_names =[j for j in ed_namesp]
    no_names = [k for k in no_namesp]
    co_names = [l for l in co_namesp]
    
    """Updates to network edges to contain all network graph information""" 
    Update_Network.update_network_edd(cod, edd, nod, ed_names, no_names, co_names)
    Update_Network.update_network_nod(cod, edd, nod, ed_names, no_names, co_names)
    print('After updating the merge events there are',len(edd), 'edges and ',
              len(nod), 'nodes in the network graph.')
    """timer"""
    if timer == True:
        print('Operation for '+str(len(cod)) +
              ' growth cones with '+ str(int(steps * instances)) +
              ' growth steps took {} min'.format((tac - tic) / 60))    
    
    """make plot"""
    if plot == True:
        Plot_Maker.plot_spline(cod, dim)
    vor = [] 
    G_vor = []
    subG_vor  = []
    """plot voronoi and dual space"""    
    if Voronoi == True:
        radius = 10
        net_obs = Network_X.voronoi_network_cluster(cod, dim, radius, 
                                                    plot_network = False,
                                                    plot_ridge_line = False) 
        
        vor, G_vor, subG_vor = net_obs
        
    G = [] 
    subG = []    
    """Build network graph via Network_X"""
    if network_x == True:
        G, subG = Network_X.create_graph_from_edd_nod(edd, nod,
                                                  plot_net = True,
                                                  plot_net_circ = False,
                                                  plot_net_spring = False)
        
    return cod, nod, edd, flavour, ed_names, no_names, co_names, vor, G_vor, subG_vor, G, subG


"""
Function name : init_manager_init()
***Description***

The function init_manager is used to create  proxy copies of nod and edd
called nodp, eddp that are manager dictionaries with the same key-value pairs
of the initial nod and edd. Also the flaovur list and the field array are copied
as proxy manager objects. flavourp is a proxy type('list') and fieldp a proxy
type('array') anly flatten type('np.ndarrays') array coulb be stored as proxy
arrays. 

***I/O***
Input parameter: 
	a) nod type('dct') dictionary that contains the initial NodeEdge class
	objects (node) 'n0','n0-1',.., 'nn','nn-1' 
	b) edd type('dict') dictionary that contains the initial NodeEdge class
	objects (edge) keys: '1',..,'n+1'
	c) flavour type('list') a list contains the range of flavour colours
	d) field type('ndarray().shape(dim)') the computation field
	
Output:
	a) mgr the manager object that creates the shared memory proxy objects
	b) eddp proxy type('dict') of edd
	c) nodp proxy type('dict') of nod
	d) flavourp proxy type('list') of flavour list
	e) fieldp proxy type('array) 1D (dim[0]*dim[1]*dim[2]) of field

Inline output:
Plot output:
Save file:
"""
def init_manager_init(edd, nod, flavour,
                      ed_names, no_names, co_names,
                      field, field_ref):
    divxp, divyp, divzp, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz= field_ref
    mgr = multiprocessing.Manager()
    eddp = mgr.dict()
    for ke,ve in edd.items():
        eddp[ke]=ve
    nodp = mgr.dict()
    for kn,vn in nod.items():
        nodp[kn]=vn
    forcep = {}
    forcep['divxp'] = divxp
    forcep['divyp'] = divyp
    forcep['divzp'] = divzp
    forcep['stpxx'] = stpxx
    forcep['stpxy'] = stpxy
    forcep['stpxz'] = stpxz
    forcep['stpyy'] = stpyy
    forcep['stpyz'] = stpyz
    forcep['stpzz'] = stpzz    
    flavourp = mgr.list(flavour)
    ed_namesp = mgr.list(ed_names)
    no_namesp = mgr.list(no_names)
    co_namesp = mgr.list(co_names)
    """mgr.Array with elements of type short interger 'h'"""	
    fieldp = mgr.Array('h',field.reshape(field.shape[0] * \
                                         field.shape[1] * \
                                         field.shape[2]))
    return mgr, eddp, nodp, flavourp, ed_namesp, no_namesp, co_namesp, fieldp, forcep

"""
Function name : init_manager_iter()
***Description***

Same function as init_manager_init without the field object.
The function was or could be used to recall loaded network objects.

---Not in use---

***I/O***
Input parameter: 
	a) nod type('dct') dictionary that contains the initial NodeEdge class
	objects (node) 'n0','n0-1',.., 'nn','nn-1' 
	b) edd type('dict') dictionary that contains the initial NodeEdge class
	objects (edge) keys: '1',..,'n+1'
	c) flavour type('list') a list contains the range of flavour colours
	
Output:
	a) mgr the manager object that creates the shared memory proxy objects
	b) eddp proxy type('dict') of edd
	c) nodp proxy type('dict') of nod
	d) flavourp proxy type('list') of flavour list

Inline output:
Plot output:
Save file:
"""
def init_manager_iter(edd, nod, flavour):
    mgr = multiprocessing.Manager()
    eddp = mgr.dict()
    for ke,ve in edd.items():
        eddp[ke]=ve
    nodp = mgr.dict()
    for kn,vn in nod.items():
        nodp[kn]=vn
    flavourp = mgr.list(flavour)
    
    return mgr, eddp, nodp, flavourp


"""
Function name : external_forces
***Description***

The function used to calculate and to call the external forces.
It works on the multilayer feature field function operators
	gradients div**
	structure tensor stp** 
Those are function operators for interpolation over the multilayer 
feature field stored in the common namespace. The function is
called in the dict_wrapper_growth function
that is the envelop for the multiproccesing function growing 
funciton of the growth process.

***I/O***
Input parameter: 
	a) cod_obj GrowthCone class object
Output:
	a) [divx, divy, divz] type('list') a list containing the three drift
        components at the GrowthCone class objects position
	b) [divxx, divxy, divxz, divyy, divyz, divzz] type('list') a list 
        containing the 6 structure tensor components at the GrowthCone class 
        objects position

Inline output:
Plot output:
Save file:
"""
def external_forces(cod_obj):
    """actual cone position"""
    x_pos,y_pos,z_pos =cod_obj.pos_list[-1]
    """
    Drift vector components
    Drift signal//all directions
    """
    divx =  forcep['divxp']([x_pos,y_pos,z_pos])[0]
    divy =  forcep['divyp']([x_pos,y_pos,z_pos])[0]
    divz =  forcep['divzp']([x_pos,y_pos,z_pos])[0]
    """
    Structure tensor components
    Diffusion metric // Diffusive tensor components
    """
    divxx = forcep['stpxx']([x_pos,y_pos,z_pos])[0]
    divxy = forcep['stpxy']([x_pos,y_pos,z_pos])[0]
    divxz = forcep['stpxz']([x_pos,y_pos,z_pos])[0]
    divyy = forcep['stpyy']([x_pos,y_pos,z_pos])[0]
    divyz = forcep['stpyz']([x_pos,y_pos,z_pos])[0]
    divzz = forcep['stpzz']([x_pos,y_pos,z_pos])[0]

    return [divx, divy, divz], [divxx, divxy, divxz, divyy, divyz, divzz]

  
"""
Function name : view_flavour_field()
***Description***

Function to generate the growth cones view in the flavour field fieldp.
The view maps the the 3x3x3 surrounding of a cone. the cone is placed
in the center.
As the flavour field, fieldp is a flatten 1D representation of field. A
coordinate transformation from 3D vectorial position np.array([x,y,z]) to
index postion [x * dim[1]*dim[2] + y *dim[2] + z]is performed. 

***I/O***
Input parameter: 
	a) cod_obj GrowthCone clasa object 
	b) dim type('ndarray().shape(3)) the dimension of the computation field 
Output:
	a) view type('np.nparray') containing references to the flavour of near
	edges. The reference is the edge's name as integer.

Inline output:
Plot output:
Save file:
"""    
def view_flavour_field(cod_obj, dim):
    
    l_pos = np.round(cod_obj.pos_list[-1]).astype('int')
    
    """The flavor field is one dimensional array -> the x,y,z pos is flatten"""
    mult_x = int((dim[1]*dim[2]))
    mult_y = int((dim[2]))
    view_axis = np.array([-2,-1,0,1,2]).astype('int')
    x_axis = (view_axis + l_pos[0]) * mult_x
    y_axis = (view_axis + l_pos[1]) * mult_y
    z_axis = (view_axis + l_pos[2])

    view= np.array([[[fieldp[min(i+j+k,int(dim[0]*dim[1]*dim[2])-1)] for i in x_axis] for j in y_axis] for k in z_axis])
      
    return view


"""
Function name : check_files_present()
***Description***

The function should be excecuted before the function init_par_growing.
It checks if all needed information and varibales are accessable to perform a
network growth simulation.
If one needed information missing an error is raised indicating what to do.

***I/O***
Input parameter: 
	a) home type('str') the name of the directory to pythrahyper_net

Output:
	a) Four path strings type('str') to the .csv- files with all information


Inline output: 'Found all parameter_files' if all neede .csv-files are present
                or an error is raised stating what to do
Plot output:
Save file:
"""
def check_files_present(home):
    
    path_internal = home + \
                    '/Parameter_Files/internal_parameters.csv'
    if Path(path_internal).is_file():
        pass
    else:
        raise Exception('You need to run the StartCondition_Maker.py ' + 
                        'in "/Start_conditions" or copy ' + 
                        'internal_parameters.csv in' + home +
                        '/pythrahyper_net/Parameter_Files/')
        
    path_growth = home + \
                  '/Parameter_Files/growth_parameters.csv'
    if Path(path_growth).is_file():
        pass
    else:
        raise Exception('You need to run the StartCondition_Maker.py ' + 
                        'in "/Start_conditions" or copy ' + 
                        'growth_parameters.csv in' + home +
                        '/pythrahyper_net/Parameter_Files/')
    
    path_startpar = home+ \
                    '/Parameter_Files/starting_positions.csv'
    if Path(path_startpar).is_file():
        pass
    else:
        raise Exception('You need to run the StartCondition_Maker.py ' +
                        'in "/Start_conditions" or copy ' + 
                        'starting_positions.csv in' + home +
                        '/pythrahyper_net/Parameter_Files/')
    
    path_multilayer = home + \
                      '/Parameter_Files/multilayer_dir_parameters.csv'
    if Path(path_multilayer).is_file():
        pass
    else:
        raise Exception('You need to run the StartCondition_Maker.py ' +
                        'and the StructField_Maker.py in "/Start_conditions"' +
                        'or copy multilayer_dir_parameter.csv in' + home +
                        '/pythrahyper_net/Parameter_Files/')
    
    print('Found all parameter_files')
    return path_multilayer, path_startpar, path_growth, path_internal  


"""
Function name : save_full_simulation()
***Description***

If this function is called information
    - GrowthCone class objects
    - NodeEdge calss objects
    - Simulation Parameter
    - multilayer computation field
    - etc
are saved. Thes save folder is created in the parent directory to 
pythrahyper_net/

***I/O***
Input parameter: 
	a) home type('str') the name of the directory to pythrahyper_net
    b) name type('str') the name of the network simulation
    c) cod type('dict') dictionary that contains all exsisting GrowthCone 
	class objects
	d) nod type('dct') dictionary that contains all exsisting NodeEdge class
	objects (node)
	e) edd type('dict') dictionary that contains all exsiting NodeEdge class
	objects (edge)
    f) fieldp proxy type('array) 1D (dim[0]*dim[1]*dim[2]) of field
    g) dim type('ndarray().shape(3)) the dimension of the computation field

Output:

Inline output:
Plot output:
Save file: the network is saved
"""
def save_full_simulation(home, network_name, cod, nod, edd, fieldp, dim):
    """
    Save flavourfield to Struct_Field_Images
    Save all network objects cod, nod, edd to save directory
    Copy all parameter files to save directory
    Copy all Struct_Field_Images *.tif to save directory
    """   
    Save_Network.save_flavourfield(fieldp, dim, home)
    """Save_directory"""
    save_path = Save_Network.create_save_directory(network_name, home)
    Save_Network.save_netwotk_cod(cod, save_path)
    Save_Network.save_netwotk_edd(edd, save_path)
    Save_Network.save_netwotk_nod(nod, save_path)
    os.popen('cp' + ' ' + path_internal + ' ' + save_path)
    os.popen('cp' + ' ' + path_growth + ' ' + save_path)
    os.popen('cp' + ' ' + path_startpar + ' ' + save_path)
    os.popen('cp' + ' ' + path_multilayer + ' ' + save_path)
    Save_Network.save_multilayer(save_path, home)



    
    
if __name__=='__main__':

    network_name = 'skynet'
    home = os.getcwd()
	
    file_paths = check_files_present(home)
    path_multilayer, path_startpar, path_growth, path_internal = file_paths
    
    """Loading and interpolation over the substrate structure, creating of 
    operator, reference functions in the global namespace"""
    field_ref = StructField_Maker.interpol_external_fields(path_multilayer)
    #divxp, divyp, divzp, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz= field_ref

    """Loading and construction of the class object dictionaries,
    global growth and envirmornt parameter"""
    obj_par_env = init_objects_parameter_enviroment(path_startpar,
                                                      path_internal,
                                                      path_growth,
                                                      path_multilayer)
    cod, nod, edd, flavour, \
    ed_names, no_names, co_names, field, steps, \
    instances, dim, radius  = obj_par_env
    
    """Constructing of the multiprocess manager and the shared memory objects
    in the global namespace"""
    mgr, eddp, nodp, flavourp, ed_namesp, \
    no_namesp, co_namesp, fieldp, forcep= init_manager_init(edd, nod, flavour,
                                                    ed_names, no_names,
                                                    co_names, field,
                                                    field_ref)
    
    """Starting the growth process """
    
    
    growing_results = init_par_growing(cod, nod, edd, steps, 
                                      instances, dim, radius,
                                      eddp, nodp, flavourp, fieldp,
                                      ed_namesp, no_namesp, co_namesp,
                                      forcep,
                                      timer = True,
                                      plot = True,
                                      Voronoi = False,
                                      network_x = True)
    
    cod, nod, edd, flavour, ed_names, no_names, co_names, \
    vor, G_vor, subG_vor, G, subG = growing_results
    
    #save_full_simulation(home, network_name, cod, nod, edd, fieldp, dim)
    #Movie_Maker.movie_maker(home, network_name, dim, cod, steps*instances)
    #Movie_Maker.spline_movie_maker(home, network_name, dim, cod, steps*instances)
   
    



         

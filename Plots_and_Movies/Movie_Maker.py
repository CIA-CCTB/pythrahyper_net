#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:20:42 2019

@author: top40ub
"""

import numpy as np
import tifffile as tif
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime
import glob
import os

"""
Function name : Movie_Maker()
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
def movie_maker(home, name, dims, cod, steps):
    homep = str(Path(home).parent)
    movie_frame = np.zeros(dims.astype('int'))

    filelist=glob.glob(homep + '/Movies/*.tif')
    for file in filelist:
        os.remove(file)
    for i in range(steps):
        for key in cod:
            if cod[key].frame_seq[0] <= i < cod[key].pos_list.shape[0]+cod[key].frame_seq[0]:
                frame_pos = cod[key].pos_list[i-cod[key].frame_seq[0]]
                frame_pos=  np.int8(frame_pos)
                movie_frame[frame_pos[0],frame_pos[1],frame_pos[2]] = 200
                
    
        tif.imsave(homep + '/Movies/' + name + str(i) + '.tif',np.uint8(movie_frame))

    return movie_frame


"""
Function name :  plot_spline_movie()
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
def plot_spline_movie(cod,frames,dim):
    plt.ioff()
    fig = plt.figure()
    ax3d = fig.add_subplot(111, projection='3d')
    i = 0
    color_map = {}
    for key in cod:
        i += 1
        r = i/len(cod)
        g = 0.5 + i/len(cod)
        if g > 1:
            g =-0.5 + i/len(cod)
        b = 1- i/len(cod)
        color_map[key]=(r,g,b)
    for key in cod.keys():
        if cod[key].frame_seq[0] <= frames:
            if cod[key].pos_list.shape[0] > 3:
                cod[key].spline_fitq()
                #ax3d.plot(cod[key].xq[:frames-cod[key].frame_seq[0]], cod[key].yq[:frames-cod[key].frame_seq[0]], cod[key].zq[:frames-cod[key].frame_seq[0]], color = color_map[key])
                
                if cod[key].sleep_merge_event == [('growing')]:
                    ax3d.plot(cod[key].new_pointsq[0][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[1][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[2][:frames-cod[key].frame_seq[0]],
                              color = 'g')
                    ax3d.scatter(cod[key].new_pointsq[0][frames-cod[key].frame_seq[0]],
                                 cod[key].new_pointsq[1][frames-cod[key].frame_seq[0]],
                                 cod[key].new_pointsq[2][frames-cod[key].frame_seq[0]],
                                 s = 1,
                                 color = 'g')
                    ax3d.scatter(cod[key].new_pointsq[0][0],
                                 cod[key].new_pointsq[1][0],
                                 cod[key].new_pointsq[2][0],
                                 s = 1,
                                 color = 'g')
                
                if frames >= cod[key].frame_seq[1]:
                    ax3d.plot(cod[key].new_pointsq[0][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[1][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[2][:frames-cod[key].frame_seq[0]],
                              color = 'r')
                    ax3d.scatter(cod[key].new_pointsq[0][-1],
                                 cod[key].new_pointsq[1][-1],
                                 cod[key].new_pointsq[2][-1],
                                 s = 1,
                                 color = 'r')
                    ax3d.scatter(cod[key].new_pointsq[0][0],
                                 cod[key].new_pointsq[1][0],
                                 cod[key].new_pointsq[2][0],
                                 s = 1,
                                 color = 'r')
                
                else:
                    ax3d.plot(cod[key].new_pointsq[0][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[1][:frames-cod[key].frame_seq[0]],
                              cod[key].new_pointsq[2][:frames-cod[key].frame_seq[0]],
                              color = 'g')
                    ax3d.scatter(cod[key].new_pointsq[0][frames-cod[key].frame_seq[0]],
                                 cod[key].new_pointsq[1][frames-cod[key].frame_seq[0]],
                                 cod[key].new_pointsq[2][frames-cod[key].frame_seq[0]],
                                 s = 1,
                                 color = 'g')
                    ax3d.scatter(cod[key].new_pointsq[0][0],
                                 cod[key].new_pointsq[1][0],
                                 cod[key].new_pointsq[2][0],
                                 s = 1,
                                 color = 'g')
                
    ax3d.set_xlim([0, dim[0]])
    ax3d.set_ylim([0, dim[1]])
    ax3d.set_zlim([0, dim[2]])
    start_pos_list = np.array([[125,125,125]])
    start_cells = list(start_pos_list)
    for i in start_cells:
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = 100 * np.outer(np.cos(u), np.sin(v)) + i[0]
        y = 100 * np.outer(np.sin(u), np.sin(v)) + i[1]
        z = 100 * np.outer(np.ones(np.size(u)), np.cos(v)) + i[2]
        ax3d.plot_surface(x, y, z, color='b', alpha= 0.1)             
    #plt.close(fig)
    return fig    

def spline_movie_maker(home, name, dims, cod, steps):
    #path_movie = create_save_directory_movie_video(name, home)
    plt.ioff()
    homep = str(Path(home).parent)
    filelist=glob.glob(homep + '/Video/*.png')
    for file in filelist:
        os.remove(file)
    for frames in range(steps):
        fig = plot_spline_movie(cod,frames,dims)
        ax =plt.gca()
        #plt.show()
        ax.grid(False)
        ax.set_axis_off()
        ax.view_init(elev=frames*0.3,azim=frames*0.6)
        plt.savefig(homep + '/Video/' + name +str(frames)+'.png', dpi=fig.dpi)
        plt.close(fig) 

"""
Function name : create_save_directory()
***Description***

The  function creates an folder in the networks/ directory with the name of the
network simultation + date_time. 
If there are no networks/ folder it is created.
The networks/ folder is placed in the parent directory home

***I/O***
Input parameter: 
	a) network_name type('str') name of the network simulation
    b) home type('str') path to the current working directory
Output:
	a) save_path type('str') path to the networks/'network_name'/ 
    folder in the parent directory of the home variable, this Networks/ folder 
    is the save folder for network simulations

Inline output:
Plot output:
Save file:
"""              
def create_save_directory_movie_video(network_name, home):
    homep = str(Path(home).parent)
    date_time = datetime.datetime.now()
    date =date_time.strftime("_%d_%b_%Y_t_%H_%M")
    
    if not os.path.exists(homep + '/Movies'):
        os.makedirs(homep + '/Movies')
    if not os.path.exists(homep + '/Movies/' + network_name):
        os.makedirs(homep + '/Movies/' + network_name + date)
    save_path_movie =homep + '/Movies/' + network_name + date
    if not os.path.exists(homep + '/Video'):
        os.makedirs(homep + '/Video')
    if not os.path.exists(homep + '/Video/' + network_name):
        os.makedirs(homep + '/Video/' + network_name + date)
    save_path_video =homep + '/Video/' + network_name + date
    return save_path_movie, save_path_video
    
if __name__=='__main__':
    movie_maker(home, network_name, dim, cod, steps*instances)
    spline_movie_maker(home, network_name, dim, cod, steps*instances)
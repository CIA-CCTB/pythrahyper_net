#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:05:12 2017

@author: top40ub
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splprep, splev
from mpl_toolkits.mplot3d import Axes3D


"""
Function name : plot_spline()
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
def plot_spline(cod,dim):
       
    
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
    for key in cod: 
        if cod[key].pos_list.shape[0] > 3:
            cod[key].spline_fitq()
            if cod[key].sleep_merge_event == [('growing')]:
                ax3d.plot(cod[key].new_pointsq[0], cod[key].new_pointsq[1], cod[key].new_pointsq[2], color = 'g', linewidth=1)
            else:
                ax3d.plot(cod[key].new_pointsq[0], cod[key].new_pointsq[1], cod[key].new_pointsq[2], color = 'r', linewidth=1)
    
    
    ax3d.set_xlim([0, dim[0]])
    ax3d.set_ylim([0, dim[1]])
    ax3d.set_zlim([0, dim[2]])
    fig.show()
        
    plt.show()
    return fig

   
if __name__=='__main__':

    plot_spline(L,dim)
         
      
        
       
        
        
        
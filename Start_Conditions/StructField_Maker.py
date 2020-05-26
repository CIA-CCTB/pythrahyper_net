#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:52:16 2019

@author: top40ub
"""
import logging
import numpy as np
import tifffile as tif
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
from time import time as time_it
import csv
from pathlib import Path
import os

'This python file loads the parameter setting for the growth simulation'
import Start_Conditions.Parameter_Importer as Parameter_Importer


"""Load computation field image, perform sturcture analysis and set the boundaries for the computation field"""
"""
Function name : load_img()
***Description***

The function takes the path string to the .csv-file containing the
path strings to the images (.tif-files) for the construction of the
multilayer computation filed.
The fields are loaded
	- the image
	- the drift image
	- the structure image
and are be default stored as float32 ndarrays.
***I/O***
Input parameter: 
	a)path_structured_images type('str') path to the .csv-file
	with the path information to the image directories and image
	names
Output:
	a) image type('ndarray').astype('float32') array to calculate
	the computation field
	b) imaged type('ndarray').astype('float32') array to calculate
	the multilayer drift field
	c) imaget type('ndarray').astype('float32') array to calculate
	the multilayer structure tensor field

Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def load_img(path_structured_images):
    
    with open(path_structured_images, 'r') as file:
        filereader = csv.reader(file)
        L = list(filereader)
        while [] in L:
            L.pop(L.index([]))
        image = np.float32(tif.imread(L[1][1]))
        imaged = np.float32(tif.imread(L[2][1]))
        imaget = np.float32(tif.imread(L[3][1]))
            
    logging.info('Loaded the images for the drift and structure layer')
    logging.info('The loaded images names are:{}'.format((L[4],L[5],L[6])) )
    return image, imaged, imaget


"""
Function name : dims_3D()
***Description***

The function takes an image array with shape(x,y,z) as imput and
returns an array containing the shape dimensions

***I/O***
Input parameter: 
	a) img tpye('ndarray').shape(x,y,z) array extract the dimensions
Output:
	a) dim type('ndarray') with shape (1,3) and the value of the 
	x,y,z dimensions

Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def dims_3D(img):
    max_x, max_y, max_z = img.shape
    dim = np.array([max_x ,max_y ,max_z])
    logging.info('The dimension of the fields are:{}'.format(([max_x, max_y ,max_z])))
    return dim


"""
Function name : drift_filter_3D()
***Description***

The function calculates 3 images for the first derivative along each Cartesiean axes (x,y,z)
of a given image type('ndarray').shape(x,y,z). These 3 images are then smooth via an isotropic
Gaussian kernel with sigma_div. The gradients are calculated with the numpy gradient operator,
the Gaussian smoothing with Gaussian fielter from scipy.ndimage.
The maximal and minimal drift value of the image are calculated.

***I/O***
Input parameter: 
	a) image type('ndarray') array to calculate the multilayer drift field
	b) sigma_div type('int') smmothing parameter for the Gaussian sigma kernel
	default for sigma is 0
Output:
	a) drift_field type('list') with 3 nd.arrays as elements, the arrays are the
	first derivatives along each axes for the imput field
	b) max_drift type('ndarray').shape(1,) with the maximal derivative value
	c) min_drift type('ndarray').shape(1,) with the mininal derivative value
	
Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def drift_filter_3D(img, sigma_div=0):
    tic = time_it()
    drift_field = np.gradient(img)
    tac = time_it()
    logging.info('The gradient calculation operation of the drift compfield took {} s to complete'.format((tac - tic) ))
    tic = time_it()
    sum1 = np.sum(drift_field[0])
    logging.info('Sum of all drift strenght of the first gradient is :{}'.format((sum1)))
    sum2 = np.sum(drift_field[1])
    logging.info('Sum of all drift strength of the second gradient is :{}'.format((sum2)))
    sum3 = np.sum(drift_field[2])
    logging.info('Sum of all drift strength of the third gradient is :{}'.format((sum3)))
    drift_field[0]= gaussian_filter(drift_field[0],sigma=sigma_div)
    sum1g = np.sum(drift_field[0])
    logging.info('Sum of all drift strength of the Gauss filtered first gradient is :{}'.format((sum1g)))
    drift_field[1]= gaussian_filter(drift_field[1],sigma=sigma_div)
    sum2g = np.sum(drift_field[1])
    logging.info('Sum of all drift strength of the Gauss filtered second gradient is :{}'.format((sum2g)))
    drift_field[2]= gaussian_filter(drift_field[2],sigma=sigma_div)
    sum3g = np.sum(drift_field[2])
    logging.info('Sum of all drift strength of the Gauss filtered first gradient is :{}'.format((sum3g)))
    tac = time_it()
    logging.info('The Gaussian filter operation and the summation of all (un-)filtered drift values for the compfield took {} s to complete'.format((tac - tic) ))
    logging.info('Sigma Gaussian filter size sigma_div:{}'.format((sigma_div)))
    tic= time_it()
    drift_strength = np.sqrt(drift_field[0]**2 +drift_field[1]**2+drift_field[2]**2)
    max_drift = np.amax(drift_strength)
    min_drift = np.amin(drift_strength)
    tac = time_it()
    logging.info('The max and min drift strength calculation operation of the compfield took {} s to complete'.format((tac - tic) ))
    logging.info('max drift:{}'.format((max_drift)))
    logging.info('min drift:{}'.format((min_drift)))
    return drift_field, max_drift, min_drift


"""
Function name : drift_filter_3D_struc_tensor()
***Description***

Same as drift_filter_3D but it is used to calcualte the derivatives for the structure
tensor.


***I/O***
Input parameter: 
	a) image type('ndarray') array to calculate the multilayer drift field
	b) sigma_div type('int') smmothing parameter for the Gaussian sigma kernel
	default for sigma is 0
Output:
	a) drift_field type('list') with 3 nd.arrays as elements, the arrays are the
	first derivatives along each axes for the imput field
	b) max_drift type('ndarray').shape(1,) with the maximal derivative value
	c) min_drift type('ndarray').shape(1,) with the mininal derivative value
	
Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def drift_filter_3D_struc_tensor(img, sigma_div=0):
    tic = time_it()
    drift_field = np.gradient(img)
    tac = time_it()
    logging.info('The gradient calculation operation for the struc tensor compfield took {} s to complete'.format((tac - tic) ))
    tic = time_it()
    sum1 = np.sum(drift_field[0])
    logging.info('Sum of all drift strenght of the first gradient is :{}'.format((sum1)))
    sum2 = np.sum(drift_field[1])
    logging.info('Sum of all drift strength of the second gradient is :{}'.format((sum2)))
    sum3 = np.sum(drift_field[2])
    logging.info('Sum of all drift strength of the third gradient is :{}'.format((sum3)))
    drift_field[0]= gaussian_filter(drift_field[0],sigma=sigma_div)
    sum1g = np.sum(drift_field[0])
    logging.info('Sum of all drift strength of the Gauss filtered first gradient is :{}'.format((sum1g)))
    drift_field[1]= gaussian_filter(drift_field[1],sigma=sigma_div)
    sum2g = np.sum(drift_field[1])
    logging.info('Sum of all drift strength of the Gauss filtered second gradient is :{}'.format((sum2g)))
    drift_field[2]= gaussian_filter(drift_field[2],sigma=sigma_div)
    sum3g = np.sum(drift_field[2])
    logging.info('Sum of all drift strength of the Gauss filtered first gradient is :{}'.format((sum3g)))
    tac = time_it()
    logging.info('The gaussian filter operation and the summation of all (un-)filtered drift values for the struc tensor compfield took {} s to complete'.format((tac - tic) ))
    logging.info('Sigma Gaussian filter:{}'.format((sigma_div)))
    tic= time_it()
    drift_strength = np.sqrt(drift_field[0]**2 +drift_field[1]**2+drift_field[2]**2)
    max_drift = np.amax(drift_strength)
    min_drift = np.amin(drift_strength)
    tac = time_it()
    logging.info('The max and min drift strength calculation operation of the struc tensor compfield took {} s to complete'.format((tac - tic) ))
    logging.info('max drift:{}'.format((max_drift)))
    logging.info('min drift:{}'.format((min_drift)))
    return drift_field, max_drift, min_drift


"""
Function name : struc_tensor_filter_3D()
***Description***


The function calculates 6 images for each of the 6 components of the structure tensor
of a given image type('ndarray').shape(x,y,z). These 6 images are then smooth via an isotropic
Gaussian kernel with sigma_div.
Basis for these six images are the 3 gradient images. Therefore these are calcualted before with
the function drift_filter_3D_struc_tensor. The 6 are the multiplication of each gradient image with
itself or with the other, the Gaussian smoothing with Gaussian fielter from scipy.ndimage.
With the nplinalg.eigvalh function the greatest and smallest eigenvalue of each possibl structure
tensor position are calculacted. 

***I/O***
Input parameter: 
	a) image type('ndarray') array to calculate the multilayer structure tensor field
	b) sigma_div2 type('int') smmothing parameter for the Gaussian sigma kernel for 
	the needed gradient images (default is 0)
	c) sigma_div2 type('int') smmothing parameter for the Gaussian sigma kernel for
	the 6 structure tensor images (default is 0)
Output:
	a) [imgxx, imgxy, imgxz, imgyy, imgyz, imgzz]  type('list') with 6 nd.arrays as elements,
	the arrays are the componetes of the structure tensor for the imput field
	b) max_eig type('ndarray').shape(1,) with the maximal eigenvalue of all structure tensors
	c) min_eig type('ndarray').shape(1,) with the minimal eigenvalue of all structure tensors
	d) eig_map type('ndarray').shape(dimx*diny*dimz,3) array with all three eigenvalues for
	each image position

Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def struc_tensor_filter_3D(img,sigma_div1=0,sigma_div2=0):
    tic = time_it()
    img, max_drift, min_drift = drift_filter_3D_struc_tensor(img,sigma_div=sigma_div1) 
    tac = time_it()
    logging.info('The gradient calculation for the struc tensor operation of the compfield took {} s to complete'.format((tac - tic) ))
    imgx = img[0]#/max_drift
    imgy = img[1]#/max_drift
    imgz = img[2]#/max_drift
    tic = time_it()
    imgxx = imgx * imgx
    imgxy = imgx * imgy
    imgxz = imgx * imgz
    imgyy = imgy * imgy
    imgyz = imgy * imgz
    imgzz = imgz * imgz
    tac = time_it()
    logging.info('The multiplication of the gradient fields took {} s to complete'.format((tac - tic) ))
    tic = time_it()
    imgxx= gaussian_filter(imgxx, sigma = sigma_div2)
    imgxy= gaussian_filter(imgxy, sigma = sigma_div2)
    imgxz= gaussian_filter(imgxz, sigma = sigma_div2)
    imgyy= gaussian_filter(imgyy, sigma = sigma_div2)
    imgyz= gaussian_filter(imgyz, sigma = sigma_div2)
    imgzz= gaussian_filter(imgzz, sigma = sigma_div2)
    tac = time_it()
    logging.info('The gauss filter operation of all combined gradient fields for the struc tensor compfield took {} s to complete'.format((tac - tic) ))
    
    eig = np.array([[imgxx.flatten(),imgxy.flatten(),imgxz.flatten()],
                    [imgxy.flatten(),imgyy.flatten(),imgyz.flatten()],
                    [imgxz.flatten(),imgyz.flatten(),imgzz.flatten()]]).T
    
    eig_map = np.linalg.eigvalsh(eig)
    """
    eig_map = np.array([np.linalg.eigvalsh([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
               for xx,xy,xz,yy,yz,zz 
               in np.nditer([imgxx,imgxy,imgxz,imgyy,imgyz,imgzz])])
    """
    logging.info('The eig_map calucaltion, calculation of the tensor distribution took {} s to complete'.format((tac - tic) ))
    tic = time_it()
    max_eig = np.amax(eig_map)
    #eig_map[eig_map <= 0.01]= 0.01
    min_eig = np.amin(eig_map)
    tac = time_it()
    logging.info('The max_eig, min_eig calculation operation of the compfield took {} s to complete'.format((tac - tic) ))
    logging.info('max_eig:{}'.format((max_eig)))
    logging.info('min_eig:{}'.format((min_eig)))
    return [imgxx, imgxy, imgxz, imgyy, imgyz, imgzz], max_eig, min_eig, eig_map


"""
Function name : inter_pol()
***Description***

The function cretes an function operator that returns an interpolated
value of a given image (img, ndarray.shape(dim[0],dim[1],dim[2]) for
any 3D position along the range for the dimension varibale (ndarray(dim[0],
dim[1),dim[2])).

***I/O***
Input parameter: 
	a) img type('ndarray').shape(dim[0],dim[1],dim[2])
	b) dim type('ndarray').shape(1,3)
Output:
	a) f_interpol the interpolation function operator for img along dim

Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def inter_pol(img,dims):

    x = np.arange(0,dims[0],1)
    y = np.arange(0,dims[1],1)
    z = np.arange(0,dims[2],1)
    f_interpol=RegularGridInterpolator((x,y,z),img,bounds_error = False, fill_value = 0)

    return f_interpol


"""
Function name : structured_field()
***Description***

The function creates the structured multilayer field images and calculates
its extrem values. The multilayered field images are stored as .tif-files and
a new .csv-file with all needed paramters, names and directories to the multilayered
field is created.

The creating starts be importing the filenames and pathes to the images used to derive
the multilayered field from a .csv-filed at postions path_structured_images.

Within in the funciton the functiuons load_img, dims_3d, drift_filter_3D,
struc_tensor_filter_3D and inter_pol are called to perform all needed calculations.
The images with the drift and structure tensor components are stored as .tif-files
and the new .csv-file in home + '/Parameter_Files/multilayer_dir_parameters.csv' is
created. the variable home is the path to the parent directory of the py.file.

***I/O***
Input parameter: 
	a) path_structured_images type('str') path to the .csv-file with the names and
	directories to the images used for calculation
	b) home type('str') directory to this .py-file
	c) sigma_divd type('int') smoothing parameter for the Gaussian kernel used
	in the creation of the drift field default for sigma is 0
	d) sigma_divt1 type('int') smoothing parameter for the Gaussian kernel used
	in the creation of the drift field for the tensor images default for sigma is 0
	d) sigma_divt2 type('int') smoothing parameter for the Gaussian kernel used
	in the creation of the tensor images field default for sigma is 0
Output:
	a) dim type('ndarray') with shape (1,3) and the value of the 
	x,y,z dimensions
	b) img1, img2 tpye('ndarray').shape(x,y,z) arrays used to calculate the multilayered
	computation field. img1 for the drift components, img2 for the structure tensor components
	c) stpxx, stpxy, stpxz, stpyy, stpyz, stpzz 6 interpolation function
	operators for the structure tensor components of the multilayered computation
	field
	d) divxp, divyp, divzp three interpoltion function operators for
	the dirft of the multilayered computation field
	e) max_drift type('ndarray').shape(1,) with the maximal derivative value
	f) min_drift type('ndarray').shape(1,) with the mininal derivative value
	g) max_eig type('ndarray').shape(1,) with the maximal eigenvalue of all structure tensors
	h) min_eig type('ndarray').shape(1,) with the minimal eigenvalue of all structure tensors
	j) eig_map type('ndarray').shape(dimx*diny*dimz,3) array with all three eigenvalues for
	each image position

Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def structured_field(path_structured_images, home, sigma_divd=0, sigma_divt1=0, sigma_divt2=0):

    logging.basicConfig(filename= home + '/Log_Files/structured_field_func.log',level=logging.DEBUG)
    logging.info('Started the structured_field function to create the multilayer computation field')
    tic=time_it()
    img, img1, img2 = load_img(path_structured_images)
    dim = dims_3D(img)
    logging.info('Starting drift_filter operation for the drift field calculation')
    img1, max_d, min_d = drift_filter_3D(img1,sigma_div=sigma_divd)
    logging.info('Finished drift_filter operation for the drift field calculation')

    logging.info('Starting struc_tensor_filter operation for the struc tensor field calculation')
    img2, max_eig, min_eig, eig_map = struc_tensor_filter_3D(img2,sigma_div1=sigma_divt1,sigma_div2=sigma_divt2)
    logging.info('Finished struc_tensor__filter operation for the struc tensor field calculation')
    tac=time_it()
    logging.info('The calculation operation took for the full 9 stack multilayer took {} s to complete'.format((tac - tic)))
    tif.imsave(home +'/Struct_Field_Images/driftfieldx.tif',np.float32(img1[0]))
    tif.imsave(home +'/Struct_Field_Images/driftfieldy.tif',np.float32(img1[1]))
    tif.imsave(home +'/Struct_Field_Images/driftfieldz.tif',np.float32(img1[2]))
    tif.imsave(home +'/Struct_Field_Images/strucxx.tif',np.float32(img2[0]))
    tif.imsave(home +'/Struct_Field_Images/strucxy.tif',np.float32(img2[1]))
    tif.imsave(home +'/Struct_Field_Images/strucxz.tif',np.float32(img2[2]))
    tif.imsave(home +'/Struct_Field_Images/strucyy.tif',np.float32(img2[3]))
    tif.imsave(home +'/Struct_Field_Images/strucyz.tif',np.float32(img2[4]))
    tif.imsave(home +'/Struct_Field_Images/struczz.tif',np.float32(img2[5]))
    logging.info('The multilayer for drift and structure the tensor parts are saved in directory ../Struct_Field_Images/ with the names driftfield*.tif and struc**.tif as float32')
    tic=time_it()
    stpxx = inter_pol(img2[0],dim)
    stpxy = inter_pol(img2[1],dim)
    stpxz = inter_pol(img2[2],dim)
    stpyy = inter_pol(img2[3],dim)
    stpyz = inter_pol(img2[4],dim)
    stpzz = inter_pol(img2[5],dim)
    divxp= inter_pol(img1[0],dim)
    divyp= inter_pol(img1[1],dim)
    divzp= inter_pol(img1[2],dim)
    tac=time_it()
    logging.info('The interpolation operation for all 9 multilayer fields took {} s to complete'.format((tac - tic) ))
    logging.info('The interpolated operator functions are called div*p and stp**')
    logging.info('Further calculated properties are the 3D dimensions dim, the min and max drift values max_d, min_d, the minimal and maximal eigenvalues of the tensors and the full distribution of the tensors eig_map')
    
    
    multilayer_dict =  {'dim':dim,
                        'max_drift':max_d,
                        'min_drift':min_d,
                        'max_eigenvalue':max_eig,
                        'min_eigenvalue':min_eig,
                        'driftx_path':home + '/Struct_Field_Images/driftfieldx.tif',
                        'drifty_path':home + '/Struct_Field_Images/driftfieldy.tif',
                        'driftz_path':home + '/Struct_Field_Images/driftfieldz.tif',
                        'tenxx_path':home + '/Struct_Field_Images/strucxx.tif',
                        'tenxy_path':home + '/Struct_Field_Images/strucxy.tif',
                        'tenxz_path':home + '/Struct_Field_Images/strucxz.tif',
                        'tenyy_path':home + '/Struct_Field_Images/strucyy.tif',
                        'tenyz_path':home + '/Struct_Field_Images/strucyz.tif',
                        'tenzz_path':home + '/Struct_Field_Images/struczz.tif'}

    path_Struct_Field_Images = home + '/Parameter_Files/multilayer_dir_parameters.csv'
    with open(path_Struct_Field_Images, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['Multilayered Field parameter and the directories of the multilayered conmputation fields'])
        for key,value in multilayer_dict.items():        
            if type(value) is np.ndarray:
                filewriter.writerow([key, *value])
            else: 
                filewriter.writerow([key, value])        


    return dim, img1, img2, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz, divxp, divyp, divzp, max_d, min_d, max_eig, min_eig, eig_map


"""
Function name : interpol_external_fields()
***Description***


The function cretes an function operator that returns an interpolated
value of a given image (img, ndarray.shape(dim[0],dim[1],dim[2]) for
any 3D position along the range for the dimension varibale (ndarray(dim[0],
dim[1),dim[2])).
The path to the images that should be interpolated over and the dimension 
variable dim are stored in a .csv-file that contains all information about
the multilayered computation field. It is called found with the path string
path_multilayer.

***I/O***
Input parameter: 
	a)path_multilayer type('str') path to the .csv-file with all information
	about the multilayered computation field
Output:
	a) divxp, divyp, divzp three interpoltion function operators for
	the dirft of the multilayered computation field
	b) stpxx, stpxy, stpxz, stpyy, stpyz, stpzz 6 interpolation function
	operators for the structure tensor components of the multilayered computation
	field


Inline output:
Plot output:
Save file: logging.info() to the logfile
"""
def interpol_external_fields(path_multilayer):
    
    multilayer = Parameter_Importer.import_multilayer_par(path_multilayer)
    dim = multilayer['dim']
    imgdx = tif.imread(multilayer['driftx_path'])
    divxp = inter_pol(imgdx,dim)
    imgdy = tif.imread(multilayer['drifty_path'])
    divyp = inter_pol(imgdy,dim)
    imgdz = tif.imread(multilayer['driftz_path'])
    divzp = inter_pol(imgdz,dim)
    imgtxx = tif.imread(multilayer['tenxx_path'])
    stpxx = inter_pol(imgtxx,dim)
    imgtxy = tif.imread(multilayer['tenxy_path'])
    stpxy = inter_pol(imgtxy,dim)
    imgtxz = tif.imread(multilayer['tenxz_path'])
    stpxz = inter_pol(imgtxz,dim)
    imgtyy = tif.imread(multilayer['tenyy_path'])
    stpyy = inter_pol(imgtyy,dim)
    imgtyz = tif.imread(multilayer['tenyz_path'])
    stpyz = inter_pol(imgtyz,dim)
    imgtzz = tif.imread(multilayer['tenzz_path'])
    stpzz = inter_pol(imgtzz,dim)
    return  divxp, divyp, divzp, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz


if __name__ == '__main__':
    home = str(Path.home())
    home = str(Path(os.getcwd()).parent)
    path_structured_images = home + '/Parameter_Files/structured_image_dir.csv'
    features = structured_field(path_structured_images, home, sigma_divd=2, sigma_divt1=2, sigma_divt2=2)
    dim, driftfield, struc_ten, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz, divxp, divyp, divzp, max_d, min_d, max_eig, min_eig, eig_map = features
    path_multilayer_par = home + '/Parameter_Files/multilayer_dir_parameters.csv'
    divxp, divyp, divzp, stpxx, stpxy, stpxz, stpyy, stpyz, stpzz= interpol_external_fields(path_multilayer_par)
    
    
    
    
    
    

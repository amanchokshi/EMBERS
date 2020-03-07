# Code developed by Jack Line and adapted here
# https://github.com/JLBLine/MWA_ORBCOMM

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import RectSphereBivariateSpline,SmoothSphereBivariateSpline

import sys
sys.path.append('../decode_rf_data')
from colormap import spectral

# Custom spectral colormap
cmap = spectral()
    
def create_model(file_name=None):
    '''Takes .ffe model, converts into healpix and smooths the response'''

    #Make an empty array for healpix projection
    beam_response = np.zeros(len_empty_healpix)*np.nan

    # Load ffe model data in, which is stored in real and imaginary components
    # of a phi/theta polaristaion representation of the beam
    # apparently this is normal to beam engineers
    data = np.loadtxt(file_name)
    theta = data[:,0]
    phi = data[:,1]
    re_theta = data[:,2]
    im_theta = data[:,3]
    re_phi = data[:,4]
    im_phi = data[:,5]

    # Convert to complex numbers
    X_theta = re_theta.astype(complex)
    X_theta.imag = im_theta

    X_phi = re_phi.astype(complex)
    X_phi.imag = im_phi


    # make a coord grid for the data in the .ffe model
    theta_range = np.linspace(epsilon,90,91)
    phi_range = np.linspace(epsilon,359.0,360)
    theta_mesh,phi_mesh = np.meshgrid(theta_range,phi_range)

    # Convert reference model into power from the complex values
    power = 10*np.log10(abs(X_theta)**2 + abs(X_phi)**2)
    
    # Make an inf values super small
    power[np.where(power == -np.inf)] = -80

    # Had a problem with edge effects, leave off near horizon values
    # Basically remove the last 90 values of power list, where θ == 90
    power = power[:-91]

    # Get things into correct shape and do an interpolation
    power.shape = phi_mesh.shape
    
    # Bivariate spline approximation over a rectangular mesh on a sphere
    # s is a paramater I had to play with to get by eye nice results
    # s: positive smoothing factor
    lut = RectSphereBivariateSpline(theta_range*(np.pi/180.0), phi_range*(np.pi/180.0), power.T,s=0.1)

    # Get the theta and phi of all healpixels
    ip = np.arange(len_empty_healpix)
    theta_rad, phi_rad = hp.pix2ang(nside, ip)

    # Evaluate the interpolated function at healpix gridpoints
    # Use spline to map beam into healpixels
    beam_response = lut.ev(theta_rad, phi_rad)

    return beam_response,theta_mesh,phi_mesh,power,theta


def plot_healpix(beam_map=None,sub=None,title=None,vmin=None,vmax=None,cmap=None):
    '''Does some plotting what what'''
    if vmin == None:
        if cmap == None:
            half_sky = hp.orthview(
                    map=beam_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,notext=True,
                    return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=beam_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,cmap=cmap,notext=True,
                    return_projected_map=True)
    else:
        if cmap == None:
            half_sky = hp.orthview(
                    map=beam_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,min=vmin,max=vmax,
                    notext=True,return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=beam_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,min=vmin,max=vmax,
                    cmap=cmap,notext=True,
                    return_projected_map=True)

    hp.graticule(dpar=10,coord='E',color='k',alpha=0.3,dmer=45)
    
    #print where(half_sky == nan)


    hp.projtext(0.0*(np.pi/180.0), 0.0, '0', coord='E')
    hp.projtext(30.0*(np.pi/180.0), 0.0, '30', coord='E')
    hp.projtext(60.0*(np.pi/180.0), 0.0, '60', coord='E')


    hp.projtext(90.0*(np.pi/180.0), 0.0, r'$0^\circ$', coord='E',color='k',verticalalignment='top')
    hp.projtext(90.0*(np.pi/180.0), 90.0*(np.pi/180.0), r'$90^\circ$', coord='E',color='k',horizontalalignment='right')
    hp.projtext(90.0*(np.pi/180.0), 180.0*(np.pi/180.0), r'$180^\circ$', coord='E',color='k')
    hp.projtext(90.0*(np.pi/180.0), 270.0*(np.pi/180.0), r'$270^\circ$', coord='E',color='k')



#Healpix constants
nside = 32
len_empty_healpix = hp.nside2npix(nside)   #12288

##Had problems with having coords = 0 when interpolating, so make a small number
##and call it epsilon for some reason

epsilon = 1.e-12

# Reference FEE models in XX & YY pols
refXX = '../../data/FEE_Reference_Models/MWA_reference_tile_FarField_XX.ffe'
refYY = '../../data/FEE_Reference_Models/MWA_reference_tile_FarField_YY.ffe'


healpix_XX,theta_mesh,phi_mesh,power_XX,theta = create_model(file_name=refXX)
healpix_YY,theta_mesh,phi_mesh,power_YY,theta = create_model(file_name=refYY)

##Plot the things to sanity check and save results
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(221,projection='polar')
ax2 = fig.add_subplot(222,projection='polar')

im1 = ax1.pcolormesh(phi_mesh*(np.pi/180.0),theta_mesh,power_XX,label='XX',vmin=-40,vmax=-20)
ax1.set_title('XX')

ax1.grid(color='k',alpha=0.3)

im2 = ax2.pcolormesh(phi_mesh*(np.pi/180.0),theta_mesh,power_YY,label='YY',vmin=-40,vmax=-20)
ax2.set_title('YY')

ax2.grid(color='k',alpha=0.3)

plot_healpix(beam_map=healpix_XX,sub=(2,2,3),title='Healpix XX',vmin=-40,vmax=-20, cmap=cmap)
plot_healpix(beam_map=healpix_YY,sub=(2,2,4),title='Healpix YY',vmin=-40,vmax=-20, cmap=cmap)

#np.savez_compressed('ref_dipole_models.npz',XX=healpix_XX, YY=healpix_YY)

#fig.savefig('reproject_dipole_models.png',bbox_inches='tight')
plt.show()

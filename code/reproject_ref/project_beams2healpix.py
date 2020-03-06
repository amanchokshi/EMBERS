# Code developed by Dr. Jack Line and adapted here
# https://github.com/JLBLine/MWA_ORBCOMM


from numpy import *
import matplotlib.pyplot as plt
import healpy as hp
from my_plotting_lib import add_colourbar
from scipy import interpolate


##Healpix constants
nside = 32
len_empty_healpix = 12288

def plot_healpix(data_map=None,sub=None,title=None,vmin=None,vmax=None,cmap=None):
    '''Does some plotting what what'''
    if vmin == None:
        if cmap == None:
            half_sky = hp.orthview(map=data_map,coord='E',half_sky=True,xsize=400,title=title,rot=(0,90,0),sub=sub,notext=True,return_projected_map=True)
        else:
            half_sky = hp.orthview(map=data_map,coord='E',half_sky=True,xsize=400,title=title,rot=(0,90,0),sub=sub,cmap=cmap,notext=True,return_projected_map=True)
    else:
        if cmap == None:
            half_sky = hp.orthview(map=data_map,coord='E',half_sky=True,xsize=400,title=title,rot=(0,90,0),sub=sub,min=vmin,max=vmax,notext=True,return_projected_map=True)
        else:
            half_sky = hp.orthview(map=data_map,coord='E',half_sky=True,xsize=400,title=title,rot=(0,90,0),sub=sub,min=vmin,max=vmax,cmap=cmap,notext=True,return_projected_map=True)
    hp.graticule(dpar=10,coord='E',color='k',alpha=0.3,dmer=45)
    #print where(half_sky == nan)


    hp.projtext(0.0*(pi/180.0), 0.0, '0', coord='E')
    hp.projtext(30.0*(pi/180.0), 0.0, '30', coord='E')
    hp.projtext(60.0*(pi/180.0), 0.0, '60', coord='E')


    hp.projtext(90.0*(pi/180.0), 0.0, r'$0^\circ$', coord='E',color='k',verticalalignment='top')
    hp.projtext(90.0*(pi/180.0), 90.0*(pi/180.0), r'$90^\circ$', coord='E',color='k',horizontalalignment='right')
    hp.projtext(90.0*(pi/180.0), 180.0*(pi/180.0), r'$180^\circ$', coord='E',color='k')
    hp.projtext(90.0*(pi/180.0), 270.0*(pi/180.0), r'$270^\circ$', coord='E',color='k')


from scipy.interpolate import RectSphereBivariateSpline,SmoothSphereBivariateSpline
##Had problems with having coords = 0 when interpolating, so make a small number
##and call it epsilon for some reason
epsilon = 1.e-12
#@profile
def create_model(file_name=None,tag=None):
    '''Takes .ffe model, converts into healpix and smooths the response'''

    ##Make an empty array for healpix projection
    beam_response = zeros(len_empty_healpix)*nan

    ##Load ffe model data in, which is stored in real and imaginary components
    ##of a phi/theta polaristaion representation of the beam
    ##apparently this is normal to beam engineers
    data = loadtxt(file_name)
    theta = data[:,0]
    phi = data[:,1]
    re_theta = data[:,2]
    im_theta = data[:,3]
    re_phi = data[:,4]
    im_phi = data[:,5]

    ##Convert to complex numbers
    X_theta = re_theta.astype(complex)
    X_theta.imag = im_theta

    X_phi = re_phi.astype(complex)
    X_phi.imag = im_phi

    ##make a coord grid for the data in the .ffe model
    theta_range = linspace(epsilon,90,91)
    phi_range = linspace(epsilon,359.0,360)
    theta_mesh,phi_mesh = meshgrid(theta_range,phi_range)

    ##Convert reference model into power from the complex values
    ##make an inf values super small
    power = 10*log10(abs(X_theta)**2 + abs(X_phi)**2)
    power[where(power == -inf)] = -80

    ##Had a problem with edge effects, leave off near horizon values
    power = power[:-91]

    ##Get things into correct shape and do an interpolation
    ##s is a paramater I had to play with to get by eye nice results
    power.shape = phi_mesh.shape

    lut = RectSphereBivariateSpline(theta_range*(pi/180.0), phi_range*(pi/180.0), power.T,s=0.1)

    ##Get the theta and phi of all healpixels
    ip = arange(len_empty_healpix)
    theta_rad, phi_rad = hp.pix2ang(nside, ip)

    ##Use spline to map beam into healpixels
    beam_response = lut.ev(theta_rad, phi_rad)

    return beam_response,theta_mesh,phi_mesh,power,theta

healpix_east,theta_mesh,phi_mesh,power_east,theta = create_model(file_name='MWA_single_dipole_finite_gnd_Eastern_FarField1.ffe',tag='Eastern')
healpix_west,theta_mesh,phi_mesh,power_west,theta = create_model(file_name='MWA_single_dipole_finite_gnd_Western_FarField1.ffe',tag='Western')

##Plot the things to sanity check and save results
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(221,projection='polar')
ax2 = fig.add_subplot(222,projection='polar')

im1 = ax1.pcolormesh(phi_mesh*(pi/180.0),theta_mesh,power_east,label='Eastern',vmin=-40,vmax=-20)
ax1.set_title('Eastern')

ax1.grid(color='k',alpha=0.3)

im2 = ax2.pcolormesh(phi_mesh*(pi/180.0),theta_mesh,power_west,label='Western',vmin=-40,vmax=-20)
ax2.set_title('Western')

ax2.grid(color='k',alpha=0.3)

plot_healpix(data_map=healpix_east,sub=(2,2,3),title='Healpix Eastern',vmin=-40,vmax=-20)
plot_healpix(data_map=healpix_west,sub=(2,2,4),title='Healpix Western',vmin=-40,vmax=-20)

savez_compressed('ung_dipole_models.npz',eastern=healpix_east,western=healpix_west)

fig.savefig('reproject_dipole_models.png',bbox_inches='tight')

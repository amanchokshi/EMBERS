# Healpix plotting code by Jack Line
# https://github.com/JLBLine/MWA_ORBCOMM

import numpy as np
import healpy as hp

def plot_healpix(data_map=None,sub=None,title=None,vmin=None,vmax=None,cmap=None):
    '''Yeesh do some healpix magic to plot the thing'''
    
    if vmin == None:
        if cmap == None:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,notext=True, return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0), sub=sub,cmap=cmap,
                    notext=True,return_projected_map=True)
    else:
        if cmap == None:
            half_sky = hp.orthview(
                    map=data_map,coord='E'
                    ,half_sky=True,xsize=400,rot=(0,90,0),
                    title=title,sub=sub,min=vmin,max=vmax,
                    notext=True,return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,rot=(0,90,0),
                    title=title,sub=sub,min=vmin,max=vmax,
                    cmap=cmap,notext=True,return_projected_map=True)

    hp.graticule(dpar=10,coord='E',color='k',alpha=0.3,dmer=45)
    
    hp.projtext(0.0*(np.pi/180.0), 0.0, '0', coord='E')
    hp.projtext(30.0*(np.pi/180.0), 0.0, '30', coord='E')
    hp.projtext(60.0*(np.pi/180.0), 0.0, '60', coord='E')


    hp.projtext(90.0*(np.pi/180.0), 0.0, r'$0^\circ$', coord='E',color='k',verticalalignment='top')
    hp.projtext(90.0*(np.pi/180.0), 90.0*(np.pi/180.0), r'$90^\circ$', coord='E',color='k',horizontalalignment='right')
    hp.projtext(90.0*(np.pi/180.0), 180.0*(np.pi/180.0), r'$180^\circ$', coord='E',color='k')
    hp.projtext(90.0*(np.pi/180.0), 270.0*(np.pi/180.0), r'$270^\circ$', coord='E',color='k')
    


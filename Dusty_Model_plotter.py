#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch, BboxConnector, BboxConnectorPatch)
plt.close('all')

# These are the most important lines of code right here. These are where you make 
# the most changes. l is your model number. it is important that it remains in 
# apostrophies as it is a string. The path is telling where your files are. That 
# will be where it goes to grab models from, and will almost certainly be different
# from your filepath


l = '0211'

path = "/home/physics/Research/DUSTY/DUSTY/Sean_mod/Model_Library/{0}/sphere-1_{0}".format(l)

# These first few functions do all of the math stuff. 
# The first one converts the magnitude into luminosity
# the second one converts the terrible units of Janskys into useful units
# the third one normalizes your function at a certain wavelength, usually at k-band
# or v-band. The important thing, 
# DO NOT TOUCH THESE   
def magnitude_conversions(L, mag, mag_0):
    f_nu = mag_0*10**(-0.4*mag)
    f_lam = f_nu*(3e8)*(1e-26)/((L**2)*0.000001)
    lam_f_lam = L*f_lam
    return f_nu, f_lam, lam_f_lam

def iso_flux_conversion(lamda, F_Jy):
    '''
    Converts units from janskys to Dusty output
    '''
    f_nu = F_Jy*10**-26
    f_lam = f_nu*300000000/(lamda*lamda*0.000001)
    return lamda*f_lam

def normalization(f_tot_model, model_flux_val, iso_flux_val):
    return f_tot_model*iso_flux_val/model_flux_val



#These three functions are only for zooming in on a section of the graph. 
#These can be removed, but it will mess up your graphs as they are 
def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = {
            **prop_lines,
            "alpha": prop_lines.get("alpha", 1) * 0.2,
            "clip_on": False,
        }

    c1 = BboxConnector(
        bbox1, bbox2, loc1=loc1a, loc2=loc2a, clip_on=False, **prop_lines)
    c2 = BboxConnector(
        bbox1, bbox2, loc1=loc1b, loc2=loc2b, clip_on=False, **prop_lines)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           clip_on=False,
                           **prop_patches)

    return c1, c2, bbox_patch1, bbox_patch2, p

def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    Connect *ax1* and *ax2*. The *xmin*-to-*xmax* range in both axes will
    be marked.

    Parameters
    ----------
    ax1
        The main axes.
    ax2
        The zoomed axes.
    xmin, xmax
        The limits of the colored area in both plot axes.
    **kwargs
        Arguments passed to the patch constructor.
    """

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, ax1.get_xaxis_transform())
    mybbox2 = TransformedBbox(bbox, ax2.get_xaxis_transform())

    prop_patches = {**kwargs, "ec": "none", "alpha": 0.2}

    c1, c2, bbox_patch1, bbox_patch2, p = connect_bbox(
        mybbox1, mybbox2,
        loc1a=3, loc2a=2, loc1b=4, loc2b=1,
        prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p

def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches = {**kwargs, "ec": "none", "alpha": 0.2}

    c1, c2, bbox_patch1, bbox_patch2, p = connect_bbox(
        mybbox1, mybbox2,
        loc1a=3, loc2a=2, loc1b=4, loc2b=1,
        prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


#this is the part you'll need to change for your specific star
#simbad magnitudes
B = 9.23
V = 7.70
G = 4.1293
J = -0.59
H = -1.55
K = -1.96

B_lam = 0.433
V_lam = 0.55
J_lam = 1.25
H_lam = 1.65
K_lam = 2.2

B_0 = 4266.7
V_0 = 3836.3
J_0 = 1594
H_0 = 1024
K_0 = 666.7




#IRAS raw, unfiltered data
#iras = 14219+2555
f12 = 8.46e02
f25 = 4.19e02
f60 = 6.92e01
f100 = 258e01


# Spectral_type = M7.5-M8
# used to constrain stellar temperature



# Spectra is now imported into your code. these files need to be in the same directory
# as your python code.
iso_wvlnth, iso_flux_Jy, iso_uncer_FJy, iso_uncer_F_norm = np.loadtxt(
    'ISO_spectra.txt', 
    unpack=True, skiprows=2)
V_date, V_mag = np.loadtxt(
    'aavsodata.txt', 
    unpack=True, usecols=(0,1), skiprows=1)
iras_wvlnth_micron, iras_lamdaF_lamda = np.loadtxt(
    'iras_spectra.txt', 
    unpack=True, skiprows=1)



# Here is where the conversions take place 
iso_flux_Wm2 = iso_flux_conversion(iso_wvlnth, iso_flux_Jy)

B_Fnu, B_flam, B_lamflam = magnitude_conversions(B_lam, B, B_0)
V_Fnu, V_flam, V_lamflam = magnitude_conversions(V_lam, V, V_0)
J_Fnu, J_flam, J_lamflam = magnitude_conversions(J_lam, J, J_0)
H_Fnu, H_flam, H_lamflam = magnitude_conversions(H_lam, H, H_0)
K_Fnu, K_flam, K_lamflam = magnitude_conversions(K_lam, K, K_0)

f12c = iso_flux_conversion(12, f12)
f25c = iso_flux_conversion(25, f25)
f60c = iso_flux_conversion(60, f60)
f100c = iso_flux_conversion(100, f100)




#Now we are loading in our dusty models. this only imports your first three models.
#They still need to be normalized

model_lamda1, f_tot_model1 = np.loadtxt(path+'.s001',  usecols=(0,1), unpack=True,)
model_lamda2, f_tot_model2 = np.loadtxt(path+'.s002',  usecols=(0,1), unpack=True,)
model_lamda3, f_tot_model3 = np.loadtxt(path+'.s003',  usecols=(0,1), unpack=True,)


# Now we are normalizing our model to a certain wavelength, in this case K-band 

norm_point = K_lamflam
model_flux_val1 = f_tot_model1[17]          #this is finding the line in the 
model_flux_val2 = f_tot_model2[17]          #dusty model that has the flux and 
model_flux_val3 = f_tot_model3[17]          #wavelength at our norm point

f_tot_norm1 = normalization(f_tot_model1, model_flux_val1, norm_point)
f_tot_norm2 = normalization(f_tot_model2, model_flux_val2, norm_point)
f_tot_norm3 = normalization(f_tot_model3, model_flux_val3, norm_point)

# Now everything is ready to plot! 

# m = 'SED Comparison - RX Boo - model {}'.format(l)
m = 'SED - RX Boo'

'''
axs = plt.figure().subplot_mosaic([
    ["zoom1"],
    ["main"],
])
axs["main"].set(xlabel='wavelength (' r'$\mu$m)', ylabel='flux',
                xscale='log',yscale='log', xlim=(0.2,200), ylim=(1e-12,1e-8))
axs["zoom1"].set(xlim=(7.5, 27.5),xscale='linear', yscale='linear', ylim=(1e-11, 4e-10), title=m)
axs["main"].plot(iras_wvlnth_micron, iras_lamdaF_lamda, 'yo', markersize=2,)
axs["main"].plot(model_lamda2, f_tot_norm2, 'k', )
axs["main"].plot(iso_wvlnth, iso_flux_Wm2, 'mo', markersize=2,)
axs["main"].plot(B_lam, B_lamflam, 'g', markersize=4, marker='.', label='BVJHK')
axs["main"].plot(V_lam, V_lamflam, 'g', markersize=4, marker='.')
axs["main"].plot(J_lam, J_lamflam, 'g', markersize=4, marker='.')
axs["main"].plot(H_lam, H_lamflam, 'g', markersize=4, marker='.')
axs["main"].plot(K_lam, K_lamflam, 'g', markersize=4, marker='.')
axs["main"].plot(12, f12c, 'r', markersize=4, marker='.', label="ISO Photometric")
axs["main"].plot(25, f25c, 'r', markersize=4, marker='.')
axs["main"].plot(60, f60c, 'r', markersize=4, marker='.')
axs["main"].plot(100, f100c, 'r', markersize=4, marker='.')
'''
fig,ax = plt.subplots()
ax.set(xscale='log', yscale='log', xlim=(1,60), ylim=(1e-11, 1e-8))
ax.set_title(m, fontsize=16)
ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
ax.set_ylabel(r'$\lambda$F$_{\lambda}$', fontsize=14)

ax.vlines(4.6, 1e-11, 1e-9, color='b')
ax.annotate('CO',
    xy=(4.6, 1e-9), xycoords='data',
    xytext=(4.0, 7e-9), textcoords='data',
    arrowprops=dict(arrowstyle="->"), color='b')

ax.vlines(6.6, 1e-11, 1e-9, color='c')
ax.annotate('H2O',
    xy=(6.6, 1e-9), xycoords='data',
    xytext=(5.3, 7e-9), textcoords='data',
    arrowprops=dict(arrowstyle="->"), color='c')

ax.vlines(8.1, 1e-11, 1e-9, color='g')
ax.annotate('SiO',
    xy=(8.1, 1e-9), xycoords='data',
    xytext=(8.1, 7e-9), textcoords='data',
    arrowprops=dict(arrowstyle="->"), color='g')

ax.plot(iras_wvlnth_micron, iras_lamdaF_lamda, 'yo', markersize=2,label='IRAS')
ax.plot(model_lamda2, f_tot_norm2, 'k', label='DUSTY')
ax.plot(iso_wvlnth, iso_flux_Wm2, 'mo', markersize=2, label='ISO')
ax.plot(B_lam, B_lamflam, 'g', markersize=4, marker='.', label='BVJHK')
ax.plot(V_lam, V_lamflam, 'g', markersize=4, marker='.')
ax.plot(J_lam, J_lamflam, 'g', markersize=4, marker='.')
ax.plot(H_lam, H_lamflam, 'g', markersize=4, marker='.')
ax.plot(K_lam, K_lamflam, 'g', markersize=4, marker='.')
ax.plot(12, f12c, 'r', markersize=4, marker='.', label="ISO Photometric")
ax.plot(25, f25c, 'r', markersize=4, marker='.')
ax.plot(60, f60c, 'r', markersize=4, marker='.')
ax.plot(100, f100c, 'r', markersize=4, marker='.')
ax.legend()
# axs["zoom1"].plot(iras_wvlnth_micron, iras_lamdaF_lamda, 'yo', markersize=2,)
# axs["zoom1"].plot(iso_wvlnth, iso_flux_Wm2, 'mo', markersize=2,)
# axs["zoom1"].plot(model_lamda2, f_tot_norm2, 'k')

# zoom_effect01(axs["zoom1"], axs["main"], 7.5, 27.5)

plt.savefig(m+'.png')
plt.show()
































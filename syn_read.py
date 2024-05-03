'''
Program to visualize the simulation data
'''
import matplotlib
matplotlib.use('TkAgg')
import Tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import MouseEvent
import numpy as np
from sepid import *
from astropy.io import fits
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from matplotlib.widgets import Cursor
from datetime import date
from datetime import datetime

# user options
#########
save = 0
verbos = 0 # print the fits headers
profile_gui = 0 # wavelength profile of clicked points
map_gui = 1 # maps of the clicked wavelength position
track = 1 # tracking the fibrils on the FOV
clear = 1 # clear canvas
g = 0.1 # gamma factor for displaying the cak image
g_integrate = 0.3 # gamma for the integrated cak map
unsharp_sigma = 0.5 # gaussian blurring radius in unsharping process
markersize = 10 # scatter plot marker size
integrate_map = 1 # display integrated map instead of the K3
zoomi, zoomn = -0.5,0.5
mean_range = 0.1 # in \AA from the line center to make the integrated map
find_k3 = 0

# Definitions
########
c = 299792458 # light speed in m/s

# DIRECTORIES
syndir_cak = '/scratch/flavio/sanja/spdeg/'
syndir_ca8 = '/storage_new/users/flavio/ca2.1024x1024x673/clv-fits/'
outdir = '/home/seki2695/OUTPUT/project3/'
datadir = '/home/seki2695/OUTPUT/project3/DATA/'

# Ca II K
cak_file = syndir_cak+file_search(syndir_cak,'ie.ca2_K_mu100.muram_plage_582800.512x512x673.fits')[0]
cak_fits = fits.open(cak_file)
cak_hdr = cak_fits[0].header
cak_flip = cak_fits[0].data # Ca II K image cube
cak = np.flip(cak_flip, axis = 0) # full-stokes ca8 in w order

nu_cak = cak_fits[1].data # Ca II K wavelength axis (in frequency)
nw = cak.shape[0] # number of wavelength positions
nx = cak.shape[1]
ny = cak.shape[2]
w_k3 = nw/2 # k3 wavelength position
w_cak = np.flip((c/nu_cak)*1e10, axis = 0) # wavelength points in Angstrum
np.save(datadir+'cak_wav', w_cak)
dw_cak = w_cak-w_cak[w_k3] # wavelength posxixtions from the line-centre
wi, wn = min(w_cak), max(w_cak) # plotting w axis range
meanf = np.max(np.where(np.round(dw_cak, decimals = 1)==mean_range
         ))
meani = np.min(np.where(np.round(dw_cak, decimals = 1)==-mean_range
         ))

# Ca 8542
ca8_file_i = syndir_ca8+file_search(syndir_ca8, 'ie.ca2_8542A_mu100.muram_plage_582800.512x512x673.fits')[0]
ca8_file_q = syndir_ca8+file_search(syndir_ca8, 'qie.ca2_8542A_mu100.muram_plage_582800.512x512x673.fits')[0]
ca8_file_u = syndir_ca8+file_search(syndir_ca8, 'uie.ca2_8542A_mu100.muram_plage_582800.512x512x673.fits')[0]
ca8_file_v = syndir_ca8+file_search(syndir_ca8, 'vie.ca2_8542A_mu100.muram_plage_582800.512x512x673.fits')[0]

ca8_fits = fits.open(ca8_file_i)
ca8_hdr = ca8_fits[0].header
ca8_i = ca8_fits[0].data # Ca 8542 image cube
ca8_q = mf.fits.getdata(ca8_file_q)
ca8_u = mf.fits.getdata(ca8_file_u)
ca8_v = mf.fits.getdata(ca8_file_v)

ca8_flip = np.stack((ca8_i, ca8_q, ca8_u, ca8_v), axis = 3) # makes the full-stokes ca8 in nu order
ca8 = np.flip(ca8_flip, axis = 0) # full-stokes ca8 in w order

nu_ca8 = ca8_fits[1].data # Ca 8542 wavelength axis (in frequency)
syn_res = cak_hdr[12] # pixel size
w_ca8 = np.flip(c/nu_ca8*1e10, axis = 0)

# Saving the cubes
if(save):
    np.save(datadir+'ca8_fullstokes_dat', ca8)
    np.save(datadir+'ca8_wav', w_ca8)
    np.save(datadir+'cak_dat', cak)
    np.save(datadir+'cak_wav', w_cak)

stop()

# print the fits header
if(verbos):
    print('Ca II K header:\n===============')
    print(repr(cak_hdr))
    print('===============')

# GUIs
plt.close("all")
f = plt.figure(figsize=(9,9), tight_layout = True) # image figure
f_profile = plt.figure(figsize=(10,4), tight_layout = True) # profile figure
ax = f.add_subplot(111) # image panel
ax_profile = f_profile.add_subplot(111) # wavelength profile panel

# wavelength axis ticks
w_ticks = np.linspace(0,nw-1,nw)
w_tick_labels = np.array(w_ticks-w_k3, dtype=int)

# x & y ticks
xticklab_no = np.array([0,10,20,30,40])
xticks = xticklab_no/syn_res
xticklabs = np.array(xticklab_no, dtype=int)

# image map
cak3 = cak[w_k3,:,:]
cak_int = np.mean(cak[meani:meanf,:,:], axis = 0)
#cak3_tile = np.tile(cak3, (3,3))
#cak_tile = np.tile(cak, (1,3,3))
#cak_int_tile = np.tile(cak_int, (3,3))

if(integrate_map==1):
    ax.imshow(gamma(unsharp(cak_int, sigma = unsharp_sigma), g_integrate), cmap = 'gray', origin = 'lower')
else:
    ax.imshow(gamma(cak3, g), cmap = 'gray', origin = 'lower')
ax.set_xlabel('x [Mm]')
ax.set_ylabel('y [Mm]')
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(10))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabs)
ax.set_yticks(xticks)
ax.set_yticklabels(xticklabs)

# wavelength profile
cak_mean = cak.mean(axis = 1).mean(axis=1) # mean profile of the whole cak map
#ax_profile.plot(dw_cak, cak_mean, color = 'gray', linestyle = '--', alpha = 0.5)
ax_profile.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
ax_profile.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax_profile.set_ylim(0,5.5e-5)
ax_profile.set_xlim(wi,wn)
#ax_profile.set_xticks(w_ticks)
#ax_profile.set_xticklabels(w_tick_labels)
ax_profile.set_ylabel('I [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$ ster$^{-1}$]')
ax_profile.set_xlabel(r'$\Delta \lambda$ [$\AA$]')
ax_profile.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)

# wavelength profile zooming panel
axins = inset_axes(ax_profile, width=4.5, height=1.3,
                   bbox_to_anchor=(0.253, 0.5),
                   bbox_transform=ax_profile.transAxes, loc=3, borderpad=0)

axins.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
axins.set_xlim(w_cak[w_k3]+zoomi,w_cak[w_k3]+zoomn)
axins.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
axins.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#stop()
if(profile_gui):
    def onclick(event):
        xx, yy = int(round(event.xdata)), int(round(event.ydata))   
        # indicating the click point
        line, = ax.plot(xx,yy, '+', markersize = markersize, markeredgewidth = 2)

        # clear canvas
        if(clear):
            # plotting the spetral profile
            ax_profile.cla() # to clear the panel
            axins.cla() # to clear the panel
            ax_profile.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            axins.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax_profile.set_ylim(0,5.5e-5)
            axins.set_xlim(w_cak[w_k3]+zoomi,w_cak[w_k3]+zoomn)
            ax_profile.set_xlim(wi,wn)
            #ax_profile.set_xticks(w_ticks)
            #ax_profile.set_xticklabels(w_tick_labels)
            ax_profile.set_ylabel('I [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$ ster$^{-1}$]')
            ax_profile.set_xlabel(r'$\Delta \lambda$ [$\AA$]')
            ax_profile.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
            axins.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
            ax_profile.plot(w_cak, cak_mean, color = 'gray', linestyle = '--', alpha = 0.25)
            axins.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
            ax_profile.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
        if (find_k3):
            if((np.argmin(cak[meani:meanf,xx, yy])!=0) and (np.argmin(cak[meani:meanf,xx, yy])!= meanf-meani)):
                k3 = meani+np.argmin(cak[meani:meanf, xx,yy])
            else:
                reng = cak[meani:meanf, xx,yy]
                for i in range (len(reng[1:-1])-3):
                    k3 = meani + np.where((reng[i+1]<reng[i+2]) and (reng[i+1]<reng[i]))
            axins.plot(w_cak[k3],cak[k3, xx%nx, yy%ny], '*', markersize = 5, markeredgewidth = 2, color = 'red')
            ax.imshow(gamma(cak[k3,:,:], g), cmap = 'gray', origin = 'lower')

        ax_profile.plot(w_cak,cak[: ,xx%nx, yy%ny], color = line.get_color(), alpha =0.25)
        axins.plot(w_cak,cak[: ,xx%nx, yy%ny], color = line.get_color(), alpha =0.25)
        f_profile.canvas.draw()
        f.canvas.draw()
    cid = f.canvas.mpl_connect('button_press_event', onclick)

    plt.show()    



if(map_gui):
    ax_profile.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
    def onclick(event):
        #dw_cak = np.array(np.linspace(0,100,101), dtype = int)
        xx, yy = event.xdata, event.ydata
        print(event.xdata, event.ydata)
        # indicating the click point
        w = np.where(np.abs(w_cak-(xx))<=np.min(np.abs(w_cak-(xx))))[0][0]
        #print(w)
        # plotting the spetral profile
        ax.cla() # to clear the panel
        cak_w = cak[w,:,:]
        ax.imshow(gamma(cak_w, g), cmap = 'gray', origin = 'lower')
        ax.set_xlabel('x [Mm]')
        ax.set_ylabel('y [Mm]')
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator(10))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabs)
        ax.set_yticks(xticks)
        ax.set_yticklabels(xticklabs)

        ax_profile.cla()
        axins.cla()
        ax_profile.plot(w_cak, cak_mean, color = 'gray', linestyle = '--', alpha = 0.5)
        ax_profile.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        axins.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax_profile.set_ylim(0,5.5e-5)
        ax_profile.set_xlim(wi,wn)
        #ax_profile.set_xticks(w_ticks)
        #ax_profile.set_xticklabels(w_tick_labels)
        ax_profile.set_ylabel('I [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$ ster$^{-1}$]')
        ax_profile.set_xlabel(r'$\lambda$ [$\AA$]')
        ax_profile.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
        axins.axvline(x = w_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
        ax_profile.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
        line, = ax_profile.plot(w_cak[w],cak_mean[w], '*', markersize = 1, markeredgewidth = 2, color = 'red')

        axins.plot(w_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')

        axins.plot(w_cak[w],cak_mean[w], '*', markersize = 1, markeredgewidth = 2, color = 'red')
        axins.set_xlim(w_cak[w_k3]+zoomi,w_cak[w_k3]+zoomn)
        f.canvas.draw()
        f_profile.canvas.draw()
        
    cid = f_profile.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

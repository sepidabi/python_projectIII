'''
Program to define bright fibril paths
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
verbos = 0 # print the fits headers
profile = 0 # plot the wavelength profile
clear = 1 # clear canvas
g = 0.1 # gamma factor for displaying the cak image
tile = 0
g_integrate = 0.2 # gamma for the integrated cak map
unsharp_sigma = 0.5 # gaussian blurring radius in unsharping process
markersize = 10 # scatter plot marker size
integrate_map = 1 # display integrated map instead of the K3
zoomi, zoomn = -0.5,0.5
mean_range = 0.1 # in \AA from the line center to make the integrated map
find_k3 = 0
#########

# DIRECTORIES
simdir_cak = '/scratch/flavio/sanja/spdeg/'
simdir_ca8 = '/storage_new/users/flavio/ca2.1024x1024x673/clv-fits/'
outdir = '/home/seki2695/OUTPUT/project3/'

# Simulation files
cak_file = simdir_cak+file_search(simdir_cak,'ie.ca2_K_mu100.muram_plage_582800.512x512x673.fits')[0]# '*ca2_K*mu100*fits')[0]
ca8_file = simdir_ca8+file_search(simdir_ca8, 'ie.ca2_8542A_mu100.muram_plage_582800.512x512x673.fits')[0]

# openning the fits file
cak_fits = fits.open(cak_file)
cak_hdr = cak_fits[0].header
cak = cak_fits[0].data # Ca II K image cube
nu_cak = cak_fits[1].data # Ca II K wavelength axis (in frequency)

ca8_fits = fits.open(ca8_file)
ca8_hdr = ca8_fits[0].header
ca8 = ca8_fits[0].data # Ca 8542 image cube
nu_ca8 = ca8_fits[1].data # Ca 8542 wavelength axis (in frequency)

sim_res = cak_hdr[12] # pixel size

nw = cak.shape[0] # number of wavelength positions
nx = cak.shape[1]
ny = cak.shape[2]
w_k3 = nw/2 # k3 wavelength position
c = 299792458 # light speed in m/s
w_cak = (c/nu_cak)*1e10 # wavelength points in Angstrum
dw_cak = w_cak-w_cak[w_k3] # wavelength positions from the line-centre
wi, wn = min(dw_cak), max(dw_cak) # plotting w axis range
meani = np.max(np.where(np.round(dw_cak, decimals = 1)==mean_range
         ))
meanf = np.min(np.where(np.round(dw_cak, decimals = 1)==-mean_range
         ))
w_ca8 = c/nu_ca8

# print the fits header
if(verbos):
    print('Ca II K header:\n===============')
    print(repr(cak_hdr))
    print('===============')

# GUIs
plt.close("all")
f = plt.figure(figsize=(9,9), tight_layout = True) # image figure
ax = f.add_subplot(111) # image panel

# wavelength axis ticks
w_ticks = np.linspace(0,nw-1,nw)
w_tick_labels = np.array(w_ticks-w_k3, dtype=int)

# x & y ticks
xticklab_no = np.array([0,10,20,30,40])
xticks = xticklab_no/sim_res
xticklabs = np.array(xticklab_no, dtype=int)

# image map
cak3 = cak[w_k3,:,:]
cak_int = np.mean(cak[meani:meanf,:,:], axis = 0)
# Tile image map
cak3_tile = np.tile(cak3, (3,3))
cak_tile = np.tile(cak, (1,3,3))
cak_int_tile = np.tile(cak_int, (3,3))

ax.set_xlabel('x [pixel]')
ax.set_ylabel('y [pixel]')
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(10))
#ax.set_xticks(xticks)
#ax.set_xticklabels(xticklabs)
#ax.set_yticks(xticks)
#ax.set_yticklabels(xticklabs)

# image display
if tile==0:
    img = cak_int
    fibril_dir = outdir+'paths/'
    fibril_fileS = file_search(fibril_dir, 'path*.npy')
else:
    img = cak_int_tile
    fibril_dir = outdir+'paths/'
    fibril_fileS = file_search(fibril_dir, 'tile_path*.npy')

ax.imshow(gamma(unsharp(img, sigma = unsharp_sigma), g_integrate), cmap = 'gray', origin = 'lower')

# overplotting the chosen fibrils
for ii in range(len(fibril_fileS)):
    fibril_file = fibril_dir+fibril_fileS[ii]
    fibril = np.load(fibril_file)
    if len(fibril)<4:
        os.system('rm '+fibril_file)
    else:
        if (np.abs(fibril[0,0]-fibril[1,0])>30 or np.abs(fibril[0,1]-fibril[1,1])>30):
            fibril = fibril[1:,:]
        ax.plot(fibril[:,0], fibril[:,1], linestyle = '--', color = 'white', linewidth = 0.75, alpha = 0.75)
        ax.text(fibril[0,0], fibril[0,1], str(ii), color = 'white')
    
if(profile):
    f_profile = plt.figure(figsize=(10,4), tight_layout = True) # profile figure
    ax_profile = f_profile.add_subplot(111) # wavelength profile panel

    # wavelength profile
    cak_mean = cak.mean(axis = 1).mean(axis=1) # mean profile of the whole cak map
    #ax_profile.plot(dw_cak, cak_mean, color = 'gray', linestyle = '--', alpha = 0.5)
    ax_profile.plot(dw_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
    ax_profile.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax_profile.set_ylim(0,5.5e-5)
    ax_profile.set_xlim(wi, wn)
    #ax_profile.set_xticks(w_ticks)
    #ax_profile.set_xticklabels(w_tick_labels)
    ax_profile.set_ylabel('I [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$ ster$^{-1}$]')
    ax_profile.set_xlabel(r'$\Delta \lambda$ [$\AA$]')
    ax_profile.axvline(x = dw_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
    
    # wavelength profile zooming panel
    axins = inset_axes(ax_profile, width=4.5, height=1.3,
                       bbox_to_anchor=(0.253, 0.5),
                       bbox_transform=ax_profile.transAxes, loc=3, borderpad=0)

    axins.plot(dw_cak, cak_mean, '.', markersize = 1., markeredgewidth = 2, color = 'gray')
    axins.set_xlim(zoomi,zoomn)
    axins.axvline(x = dw_cak[w_k3], ymin=0,ymax = 5,linestyle = ':', color = 'gray', alpha = 0.5)
    axins.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


coord = []

# GUI stuff
cursor = Cursor(ax, horizOn=True, vertOn=True, useblit=True,
                color = 'r', linewidth = 0.75, alpha = 0.75)
    
def onclick(event):
    xx, yy = int(round(event.xdata))%nx, int(round(event.ydata))%ny
    global coord
    coord.append((xx, yy))
    ax.plot(xx,yy, '+', color = 'red')
    print('x=%1.2f, y=%1.2f' %
          (np.round(xx, decimals = 2),
           np.round(yy, decimals = 2)))
    f.canvas.draw() #redraw the figure

c = 1
while c==1:
    zoom_ok = False
    print('\n1. Zoom or pan to view, \n2. Press spacebar and click along the path.\n ')
    while not zoom_ok:
        zoom_ok = plt.waitforbuttonpress()
        
    coord = []
    cid = f.canvas.mpl_connect('button_press_event', onclick)
    plt.show()    

    print('\n3. press "c" When done with one path.')
    stop()
    
    # save the clicked points
    now = datetime.now()
    #osc_fname = outdir+"/paths/path"+now.strftime("%y%m%d-%H%M%S")+".txt"        
    #np.savetxt(osc_fname, coord, fmt='%3.8f', delimiter=' ', newline='\n', header='x [pixel], y [pixel]', footer='', comments='# ', encoding=None)
    if tile==0:
        fname = outdir+"/paths/path"+now.strftime("%y%m%d-%H%M%S")+".npy"
    else:
        fname = outdir+"/paths/tile_path"+now.strftime("%y%m%d-%H%M%S")+".npy"
        
    np.save(fname, np.array(coord))
        
    c = input("Still tracking? (0 or 1): ")
    
    if c==0:
        break



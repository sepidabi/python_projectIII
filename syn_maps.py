import sparsetools as sp
import numpy as np
import matplotlib.pyplot as plt
import m3d
import witt
import tqdm
import time
from sepid import *
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# user options
#########
save = 0
verbos = 0 # print the fits headers
profile_gui = 1 # wavelength profile of clicked points
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
res = 25.575

# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
fontsize = 8.

# FUNCTIONS
# Conversionof intensity in nu to lambda
def Iwav(Inu, nu, c = 299792458): 
   Iwav = np.zeros(Inu.shape)
   for i in range(Inu.shape[1]):
      for j in range(Inu.shape[2]):
         Iwav[:,i,j] = (nu**2/c)*Inu[:,i,j]
   return Iwav

# Conversion of nu to lambda
def w(nu, c = 299792458):
   w = c/nu
   return w

# DIRECTORIES
syndir_cak = '/scratch/flavio/sanja/spdeg/'
syndir_ca8 = '/storage_new/users/flavio/ca2.1024x1024x673/clv-fits/'
outdir = '/home/seki2695/OUTPUT/project3/'
datadir = '/home/seki2695/OUTPUT/project3/DATA/'

# Ca II K
cak_file = syndir_cak+file_search(syndir_cak,'ie.ca2_K_mu100.muram_plage_582800.512x512x673.fits')[0]
cak_fits = fits.open(cak_file)
cak_hdr = cak_fits[0].header
nu_cak = cak_fits[1].data # Ca II K wavelength axis (in frequency)
cak_flip = Iwav(cak_fits[0].data, nu_cak)  # Ca II K image cube
cak = np.flip(cak_flip, axis = 0) # full-stokes ca8 in w order

nw = cak.shape[0] # number of wavelength positions
nx = cak.shape[1]
ny = cak.shape[2]
w_k3 = nw/2 # k3 wavelength position
w_cak = np.flip((c/nu_cak)*1e10, axis = 0) # wavelength points in Angstrum
np.save(datadir+'cak_wav', w_cak)
dw_cak = w_cak-w_cak[w_k3] # wavelength positions from the line-centre
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
nu_ca8 = ca8_fits[1].data # Ca 8542 wavelength axis (in frequency)
ca8_i = ca8_fits[0].data # Ca 8542 image cube
ca8_q = mf.fits.getdata(ca8_file_q)
ca8_u = mf.fits.getdata(ca8_file_u)
ca8_v = mf.fits.getdata(ca8_file_v)

ca8_flip = np.stack((Iwav(ca8_i, nu_ca8), ca8_q, ca8_u, ca8_v), axis = 3) # makes the full-stokes ca8 in nu order
ca8 = np.flip(ca8_flip, axis = 0) # full-stokes ca8 in w order

syn_res = cak_hdr[12] # pixel size
w_ca8 = np.flip(c/nu_ca8, axis = 0)

# Fe I 6301 & 2
fe_dir = '/home/seki2695/INV/stic/III/stic_new/synth/'
fe_file = sp.profile(fe_dir+'synthetic_taucrop_580000_pg_feI.nc')
fe = fe_file.dat[0,:,:,:,:]
fe_cont = fe[:,:,-1,0]
w_fe = fe_file.wav

# photosphere LPtot
LPtot = np.zeros((fe.shape[0], fe.shape[1]))
CPtot = np.zeros((fe.shape[0], fe.shape[1]))
CPt_l = np.zeros((fe.shape[0], fe.shape[1]))
CPt_r = np.zeros((fe.shape[0], fe.shape[1]))
StV_l = np.zeros((fe.shape[0], fe.shape[1]))
StV_r = np.zeros((fe.shape[0], fe.shape[1]))
for ii in range(fe.shape[2]-1):
    LPtot += np.sqrt(fe[:,:,ii,1]**2 + fe[:,:,ii,2]**2)#/fe_cont
LPtot = LPtot/fe.shape[2]
# photosphere CPtot
for l in range(42,49):
   delt = 1
   StV_l += fe[:,:,l,3]
CPt_l = CPt_l + delt*StV_l
for r in range(54,65):
   delt = -1
   StV_r += fe[:,:,r,3]
CPt_r = CPt_r + delt*StV_r

CPtot = (CPt_l+CPt_r)#/fe_cont

# chromosphere LPtot
LPtotc = np.zeros((fe.shape[0], fe.shape[1]))
CPtotc = np.zeros((fe.shape[0], fe.shape[1]))
CPt_lc = np.zeros((fe.shape[0], fe.shape[1]))
CPt_rc = np.zeros((fe.shape[0], fe.shape[1]))
StV_lc = np.zeros((fe.shape[0], fe.shape[1]))
StV_rc = np.zeros((fe.shape[0], fe.shape[1]))
for ii in range(ca8.shape[0]):
    LPtotc += np.sqrt(ca8[ii,:,:,1]**2 + ca8[ii,:,:,2]**2)#/fe_cont
LPtotc = LPtotc/ca8.shape[0]
# chromospere CPtot
for l in range(9,21):
   delt = 1
   StV_lc += ca8_v[l,:,:]
CPt_lc = CPt_lc + delt*StV_lc

for r in range(31,38):
   delt = -1
   StV_rc += ca8_v[r,:,:]
CPt_rc = CPt_lc + delt*StV_lc

CPtotc = (CPt_lc+CPt_rc)


#stop()

# continuum at 4000 AA
cont_dir = fe_dir
cont_file = sp.profile(cont_dir+'cont_4000.nc')
cont = cont_file.dat[0,:,:,0,0]
w_cont = cont_file.wav[0]

#stop()

# Plotting the maps
plt.close('all')
f , ax= plt.subplots(ncols = 3, nrows =3, figsize = (8,8))

# x & y ticks
xticklab_no = np.array([0,10,20,30])
xticks = xticklab_no*res
xticklabs = np.array(xticklab_no, dtype=int)

maps = np.array(([cont, gamma(LPtot, g = 0.35), CPtot*1.e5],
                 [gamma(ca8[51/2, :,:,0], g = 0.4), gamma(LPtotc, g = 0.44), CPtotc],
                 [gamma(cak[101/2+20,:,:], g = 0.15), gamma(cak[101/2,:,:,], g = 0.15), gamma(cak[101/2-20,:,:], g = 0.15)]))

cmap = [['gray', 'gray', 'gray'],
        ['gray', 'gray', 'gray'],
        ['gray', 'gray', 'gray']]

vmin = np.array([[np.min(maps[0,0,:,:]), np.min(maps[0,1,:,:]), -10.],
                 [np.min(maps[1,0,:,:]), np.min(maps[1,1,:,:]), -0.5],
                 [np.min(maps[2,0,:,:]), np.min(maps[2,1,:,:]), np.min(maps[2,2,:,:])]])

vmax = np.array([[np.max(maps[0,0,:,:]), np.max(maps[0,1,:,:]), 10.],
                 [np.max(maps[1,0,:,:]), np.max(maps[1,1,:,:]), 0.5],
                 [np.max(maps[2,0,:,:]), np.max(maps[2,1,:,:]), np.max(maps[2,2,:,:])]])

title = [[r'a) continuum 4000 $\mathrm{\AA}$',r'b) Fe I 6301 LP', r'c) Fe I 6301 CP'],
[r'd) Ca II 8542 $\mathrm{\AA}$ core',r'e) Ca II 8542 LP', r'f) Ca II 8542 CP'],
[r'g) C II K core',r'h) Ca II K$_{\mathrm{2V}}$', r'i) Ca II K$_\mathrm{{2R}}$']]

for i in range(3):
   for j in range(3):
      map = maps[i,j,:,:]
      factor = 2.5
      img = ax[i,j].imshow(map, origin = 'lower', cmap = cmap[i][j], vmin = vmin[i,j]
                           , vmax = vmax[i,j]
                           )
      '''axin = inset_axes(ax[i,j], #axis
                             width="3%",  # width = 10% of parent_bbox width
                             height="90%",  # height : 50%
                             loc='center left',
                             bbox_to_anchor=(1.02, 0., 1, 1),
                             bbox_transform=ax[i,j].transAxes,
                             borderpad=0,
      )

      cb = plt.colorbar(img, cax = axin, orientation="vertical")
      cb.ax.tick_params(labelsize=8) 
      cb.set_label(title[i][j], rotation=270, labelpad=15)
      '''
      ax[i,j].set_xticks(xticks)
      ax[i,j].set_yticks(xticks)
      ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(10))
      ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(10))
      ax[i,j].tick_params(axis="both",which='both',direction="in",labelsize = fontsize, color = 'white')
      ax[i,j].spines['bottom'].set_color('white')
      ax[i,j].spines['top'].set_color('white') 
      ax[i,j].spines['right'].set_color('white')
      ax[i,j].spines['left'].set_color('white')
      ax[i,j].text(30,30, title[i][j], color = 'white')

      if j==0:
         ax[i,j].set_ylabel('y [Mm]', fontdict = font)
         ax[i,j].set_yticklabels(xticklabs)
      else:
         ax[i,j].set_yticklabels([])

      if (i==2 and j<=2):
         ax[i,j].set_xlabel('x [Mm]', fontdict = font)
         ax[i,j].set_xticklabels(xticklabs)
      else:
         ax[i,j].set_xticklabels([])

plt.subplots_adjust(left=0.055
                                  , bottom=0.05
                                  , right=0.99
                                  , top=0.99
                                  , wspace=0.0
                                  , hspace=0.0
)

plt.show()

filename = outdir + 'syn_maps.pdf'
f.savefig(filename, dpi = 1000)
print 'file saved to: '+ filename

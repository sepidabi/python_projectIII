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

# user definitions
tau_scl=1

#Definitions
snapnumber='580000'
dir='/home/seki2695/INV/stic/III/sanja/'
outdir = '/home/seki2695/OUTPUT/project3/'
s = (1024,1024,673)
res = 25.575
zres = 21.4844/672
zn1 = 648 # logt ~ 0
zn2 = 605 # logt ~ -4.3
nx=s[0]; ny=s[1]; ndep=s[2]; nt=1
gs=nx*ny*ndep*4
# fonts
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8.,
        }
fontsize = 8.

# plasma parameters
fname='atm3d.muram_plage_580000.1024x1024x673'
temp = np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*1, order='F').transpose((1,0,2))
vlos  = np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*4, order='F').transpose((1,0,2))
rho = np.memmap(dir+fname, dtype='float32', mode='r',shape=s, offset=gs*5, order='F').transpose((1,0,2))

# tau info
fname = "%s/tau500"  % (dir)
tau=np.memmap(fname, dtype='float32', mode='r',shape=s, order='F').transpose((1,0,2))
ltau = np.log10(tau)
# electron density
fname = '%s/ne3d-LTE.muram_plage_580000.1024x1024x673' % (dir)
nne = np.memmap(fname,dtype='float32', mode='r',shape=s,order='F').transpose((1,0,2))

# magnetic field
fname='magnetic.muram_plage_580000.1024x1024x673'
bxx=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*0, order='F').transpose((1,0,2))
byy=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*1, order='F').transpose((1,0,2))
Bln=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*2, order='F').transpose((1,0,2))


Bho = (bxx**2+byy**2)**0.5
azi = np.arctan2(np.abs(byy), bxx)
Btot = np.sqrt(Bho**2+Bln**2)

# height in the atmosphere
height = np.zeros((nx,ny,ndep))
z = np.arange(0,1024, dtype=np.float)-377.
#print(z.shape) 
z=z[::-1]*21.4844*1.e5
z=z[315:675]
#print(z.size)
for ii in range(360):
   height[:,:,ii] = z[ii]


if (tau_scl==1):

   temp_0 = np.zeros((nx,ny))
   vlos_0 = np.zeros((nx,ny))
   Bln_0 = np.zeros((nx,ny))
   temp_n = np.zeros((nx,ny))
   vlos_n = np.zeros((nx,ny))
   Bln_n = np.zeros((nx,ny))

   t0 = 0.
   tn = -4.31

   avg_fact = 6
   
   for i in range(nx):
      for j in range(ny):
         tau_0 = np.where(np.abs(ltau[i,j,:]-t0)==np.min(np.abs(ltau[i,j,:]-t0)))[0][0]
         temp_0[i,j] = np.mean(temp[i,j,tau_0-avg_fact:tau_0+avg_fact])
         vlos_0[i,j] = np.mean(vlos[i,j,tau_0-avg_fact:tau_0+avg_fact])
         Bln_0[i,j] = np.mean(Bln[i,j,tau_0-avg_fact:tau_0+avg_fact])

         tau_n = np.where(np.abs(ltau[i,j,:]-tn)==np.min(np.abs(ltau[i,j,:]-tn)))[0][0]
         temp_n[i,j] = np.mean(temp[i,j,tau_n-avg_fact:tau_n+avg_fact])
         vlos_n[i,j] = np.mean(vlos[i,j,tau_n-avg_fact:tau_n+avg_fact])
         Bln_n[i,j] = np.mean(Bln[i,j,tau_n-avg_fact:tau_n+avg_fact])
#print([temp_0[i,j], vlos_0[i,j], Bln_0[i,j]])

   maps = np.array(([temp_0*1e-3, temp_n*1e-3],
                    [-vlos_0, -vlos_n],
                    [Bln_0*1e-3, Bln_n*1e-3]))
   
   cmap = [['gist_heat', 'gist_heat'],
           ['bwr', 'bwr'],
           ['RdGy_r', 'RdGy_r']]
   
   vmin = np.array(([5.2, 4.],
                    [-5.,-9.5],
                    [-2,-1]))
   
   vmax = np.array(([7.9, 8],
                    [5., 9.5],
                    [2.,1.]))
   
   title = [['',r'$T$ [kK]'],
            ['',r'$v_{\mathrm{LOS}}$ [km s$^{-1}$]'],
            ['',r'$B_{\perp}$ [kG]']]

   tit_0 = r'log($\tau_{\mathrm{500}}$) = 0.0'
   tit_n = r'log($\tau_{\mathrm{500}}$) = -4.3'

   filename = outdir + 'sim_maps_ltau.pdf'

if(tau_scl==0):
   maps = np.array(([temp[:,:,zn1]*1e-3, temp[:,:,zn2]*1e-3],
                    [-vlos[:,:,zn1], -vlos[:,:,zn2]],
                    [Bln[:,:,zn1]*1e-3, Bln[:,:,zn2]*1e-3]))
   
   cmap = [['hot', 'hot'],
           ['bwr', 'bwr'],
           ['RdGy_r', 'RdGy_r']]
   
   vmin = np.array(([4., 2.],
                    [-5.5,-10],
                    [-2,-1]))
   
   vmax = np.array(([11., 10],
                    [5.5, 10],
                    [2.,1.]))
   tit_0 = r'z = '+ str(np.round(21.4844-zn1*zres, decimals = 2)) +' Mm'
   tit_n = r'z = '+ str(np.round(21.4844-zn2*zres, decimals = 2)) +' Mm'

   filename = outdir + 'sim_maps_z.pdf'

   
# Plotting the maps
plt.close('all')
f , ax= plt.subplots(ncols = 2, nrows =3, figsize = (8.,10.))

f.subplots_adjust(left=0.02
                                  , bottom=0.04
                                  , right=0.935
                                  , top=0.97
                                  , wspace=0.01
                                  , hspace=0.03
)

# x & y ticks
xticklab_no = np.array([0,10,20,30,40])
xticks = xticklab_no*res
xticklabs = np.array(xticklab_no, dtype=int)
   
for i in range(3):
   for j in range(2):
      map = maps[i,j,:,:]
      factor = 2.5
      img = ax[i,j].imshow(map, origin = 'lower', cmap = cmap[i][j], vmin = vmin[i,j], vmax = vmax[i,j], interpolation = 'hermite')
      axin = inset_axes(ax[i,j], #axis
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
      
      ax[i,j].set_xticks(xticks)
      ax[i,j].set_yticks(xticks)
      ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(10))
      ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(10))
      ax[i,j].tick_params(labelsize = fontsize)

      if j==0:
         ax[i,j].set_ylabel('y [Mm]', fontdict = font)
         ax[i,j].set_yticklabels(xticklabs)
      else:
         ax[i,j].set_yticklabels([])

      if (i==2 and j<=1):
         ax[i,j].set_xlabel('x [Mm]', fontdict = font)
         ax[i,j].set_xticklabels(xticklabs)
      else:
         ax[i,j].set_xticklabels([])

      if (i==0 and j==0):
         ax[i,j].set_title(tit_0)
      if (i==0 and j==1):
         ax[i,j].set_title(tit_n)
         

plt.show()

f.savefig(filename, dpi = 1000)
print 'file saved to: '+ filename


''' This script is to prepare cuts from the synthesied maps
to perform inversions on'''

from sepid import *
import mfits as mf
import sparsetools as sp
import matplotlib.patches as patches

# definitions
syn_dir = '/home/seki2695/INV/stic/III/calibrated_cubes/'

# cubes
fe = mf.readfits(syn_dir+'MuRAMsyn_fe_6302_calibrated.fits')
ca8 = mf.readfits(syn_dir+'MuRAMsyn_ca8_8542_calibrated.fits')
cak = mf.readfits(syn_dir+'MuRAMsyn_cak_3950_calibrated.fits')

# getting the cut coords
plt.close("all")
fig, ax = plt.subplots(1,3,figsize = (24,8))
w = [5,-6] # wavelength integration range
ax[0].plot(cak[0,:-1,0,0], np.mean(np.mean(cak[1,:-1,:,:], axis = 1), axis = 1), linestyle = '--', color = 'black')
ax[0].plot(cak[0,:-1,0,0], np.mean(np.mean(cak[1,:-1,:,:], axis = 1), axis = 1), 'r.')
ax[1].imshow(gamma(np.mean(cak[1,w[0]:w[1],:,:], axis = 0), g = 0.1), origin = 'lower', cmap = 'gray')
ax[2].imshow(cak[1,-1,:,:,], cmap = 'gray', origin = 'lower')

plt.show()
plt.tight_layout()

x0, dx = input('x0 = '),input('dx = ')
y0, dy = input('y0 = '), input('dy = ')
xn, yn = x0+dx, y0+dy

box = patches.Rectangle((x0,y0), dx, dy, linewidth=1, edgecolor='r', facecolor='none')
ax[1].add_patch(box)
#ax[2].add_patch(boxx)

plt.show()
stop()
mf.writefits(syn_dir+'fe_test1.fits', fe[:,:,x0:xn, y0:yn])
mf.writefits(syn_dir+'ca8_test1.fits', ca8[:,:,x0:xn, y0:yn])
mf.writefits(syn_dir+'cak_test1.fits', cak[:,:,x0:xn, y0:yn])

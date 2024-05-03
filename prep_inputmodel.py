import sparsetools as sp
import numpy as np
import matplotlib.pyplot as plt
import m3d
import witt
import tqdm
import time
from sepid import *

# tau range from previous inversions
#taumin = -7.8
#taumax= 1.0
#dtau = 0.15 # changed from 0.2, to match the 1st cycle with finer grid
#ntau = int((taumax-taumin)/dtau) + 1
#tau = np.arange(ntau, dtype='float64')/(ntau-1.0) * (taumax-taumin) + taumin

#z_new = z[330:780]
snapnumber='580000'
#dir='/srv/obelix_3/sdani/new_runs/new_fem_hires/synth/'+snapnumber+'nlte/'
dir='/home/seki2695/INV/stic/III/sanja/'
s = (1024,1024,673)
#s_new = (1024,1024,ntau)

nx=s[0]; ny=s[1]; ndep=s[2]; nt=1
#ndep_new = ntau
astic = sp.model(nx=nx, ny=ny, ndep=ndep, nt=nt)
#astic_new = sp.model(nx=nx, ny=ny, ndep=ndep_new, nt=nt)
gs=nx*ny*ndep*4
#gs_new=nx*ny*ntau*4

# plasma parameters
fname='atm3d.muram_plage_580000.1024x1024x673'
astic.temp[0,:,:,:]  = np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*1, order='F').transpose((1,0,2))
astic.vlos[0,:,:,:]  = np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*4, order='F').transpose((1,0,2))
astic.rho[0,:,:,:]   = np.memmap(dir+fname, dtype='float32', mode='r',shape=s, offset=gs*5, order='F').transpose((1,0,2))


# tau info
fname = "%s/tau500"  % (dir)
bla=np.memmap(fname, dtype='float32', mode='r',shape=s, order='F').transpose((1,0,2))
astic.ltau[0,:,:,:]=np.log10(bla)

# electron density
fname = '%s/ne3d-LTE.muram_plage_580000.1024x1024x673' % (dir)
bla = np.memmap(fname,dtype='float32', mode='r',shape=s,order='F').transpose((1,0,2))
astic.nne[0,:,:,:]=bla


# magnetic field
fname='magnetic.muram_plage_580000.1024x1024x673'
bxx=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*0, order='F').transpose((1,0,2))
bzz=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*1, order='F').transpose((1,0,2))
astic.Bln[0,:,:,:]=np.memmap(dir+fname, dtype='float32', mode='r',shape=s,  offset=gs*2, order='F').transpose((1,0,2))


astic.Bho[0,:,:,:] = (bxx**2+bzz**2)**0.5
astic.azi[0,:,:,:] = np.arctan2(np.abs(bzz), bxx)

# height in the atmosphere
z = np.arange(0,1024, dtype=np.float)-377.
print(z.shape) 
z=z[::-1]*21.4844*1.e5
z=z[315:675]
print(z.size)
for ii in range(360):
   astic.z[:,:,:,ii] = z[ii]

#stop()

witt = witt.witt()
for i in range(nx):
   for j in range(ny):
      for t in range(ndep):
         astic.pgas[0,i,j,t] = witt.pg_from_rho(astic.temp[0,i,j,t], astic.rho[0,i,j,t])
   print(str(i)+'/'+str(nx))


# Cropping the tau and interpolating it to the old observations tau grid
#for i in range(nx-1):
 #  for j in range(ny-1):
  #    astic_new.ltau[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.ltau[0,i,j,:])
   #   astic_new.temp[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.temp[0,i,j,:])
    #  astic_new.vlos[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.vlos[0,i,j,:])
     # astic_new.rho[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.rho[0,i,j,:])
     # astic_new.nne[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.nne[0,i,j,:])
     # astic_new.Bln[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.Bln[0,i,j,:])
     # astic_new.Bho[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.Bho[0,i,j,:])
     # astic_new.azi[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.azi[0,i,j,:])
     # astic_new.z[0,i,j,:] = np.interp(tau, astic.ltau[0,0,0,:], astic.z[0,i,j,:])

        
    #a1 = a.extract(x0=1050, x1=1600, z0=680, z1=1423)
#a1 = a.extract(x0=320, x1=670, z0=0, z1=289, y0=0, y1=511)
#a1 = astic.extract(x0=411, x1=511, z0=0, z1=360, y0=0, y1=100)
#a1 = astic.extract(x0=200, x1=500, z0=0, z1=360, y0=150, y1=400)
#a1 = m.extract(x0=250, x1=350, z0=0, z1=360, y0=300, y1=400)
#a1.write(dir2+'modelin_plage_tau_595000_cropped.nc' , write_all=True)
#a1.cmass[:,:,:,:] = S.get_cmass(a1.z, a1.rho, a1.pgas, a1.temp)
    #a1 = S.optimize_dep(a1)


astic.write(dir+'modelin_tau_'+snapnumber+'.nc' , write_all=True) # full tau grid
#astic_new.write(dir+'modelin_new_tau_'+snapnumber+'.nc' , write_all=True) # inversion tau grid
    
    
print("test3")

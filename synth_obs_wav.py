'''This script is to resample the synthesized Fe I 6301&2
to the wavelength sampling of the observations of paper I'''

from sepid import *
import mfits as mf
import sparsetools as sp

# user option
resample = 0

# definitions
syn_dir = '/home/seki2695/INV/stic/III/stic_new/synth/'
synout_dir = '/home/seki2695/OUTPUT/project3/DATA/'
syn_file = sp.profile(syn_dir + 'synthetic_taucrop_580000_pg_feI.nc')
obs_dir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
obs_file = mf.readfits(obs_dir+'b_obs6302_06_165307.fits')

syn_dat = syn_file.dat[:,:,:,1:,:]
syn_wav = syn_file.wav[1:]
stop()
obs_wav = obs_file[0,:,0]
obs_dat = obs_file[1:,:,:]
nw_obs = len(obs_wav)

nt = syn_dat.shape[0]
nx = syn_dat.shape[1]
ny = syn_dat.shape[2]

w = [0,-1]#[np.where(np.round(syn_wav, decimals = 2) == np.round(obs_wav, decimals = 2)[0])[0][0],
     #np.where(np.round(syn_wav, decimals = 2) == np.round(obs_wav, decimals = 2)[-1])[0][0]]

syn_obs_dat = np.zeros((len(syn_wav[w[0]:w[1]]), nx,ny,4))#np.zeros((nt, nx, ny, nw_obs, 4))
syn_obs_wav = syn_wav[w[0]:w[1]]

if (resample==1):
    # resampling to observations
    for i in range(nw_obs):
        idx = np.where(np.round(syn_wav, decimals = 2) == np.round(obs_wav, decimals = 2)[i])[0][0]
        syn_obs_dat[i,:,:,:] = syn_dat[0,:,:,idx,:]
        filename = synout_dir+'fe_fullstokes_dat_resampled.npy'
        filename_w = synout_dir+'fe_wav_resampled.npy'

else:
    # making compatible npy cube
    for i in range(len(syn_wav[w[0]:w[1]])):
        syn_obs_dat[i,:,:,:] = syn_dat[0,:,:,w[0]+i,:]
    filename = synout_dir+'fe_fullstokes_dat.npy'
    filename_w = synout_dir+'fe_wav.npy'
        
filename_cont = synout_dir+'cont_4000.11.npy'
np.save(filename, syn_obs_dat)
np.save(filename_w, syn_obs_wav)
np.save(filename_cont, syn_file.dat[0,:,:,0,0])
        

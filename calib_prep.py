'''This script is to calibrate the intensity
of the synthesised profiles'''

from sepid import *
import mfits as mf
import sparsetools as sp
import ISPy.spec.calib as calib
plt.close("all")

calibrate = 1
interpolate = 0 # interpolate to obseravations wavelength point

# User options
line = 6302 # 6302/8542/3950
save_fits = 0
x = [0,285]
y = [477,778]

# Definitions
syn_dir = '/home/seki2695/OUTPUT/project3/DATA/'
obs_dir = '/scratch/sepid/DATA/AR/plage/2016.09.15/OUTPUT/'
c = 299792458 # light speed in m/s

# Fe 6301 & 2
if (line==6302):

    # definitions
    obs_file = mf.readfits(obs_dir+'b_obs6302_06_165307.fits')
    
    syn_dat = np.load(syn_dir + 'fe_fullstokes_dat.npy')
    syn_wav = np.load(syn_dir + 'fe_wav.npy')
    fe_i = syn_dat[:,:,:,0]
    fe_q = syn_dat[:,:,:,1]
    fe_u = syn_dat[:,:,:,2]
    fe_v = syn_dat[:,:,:,3]
    
    obs_wav = syn_wav
    obs_dat = obs_file[0:,:,:]
    obs_dat_avg = np.mean(obs_dat[1,:,:], axis = 1)
    nw_obs = len(obs_wav)
    
    nt = 1
    nx = syn_dat.shape[1]
    ny = syn_dat.shape[2]
    
    syn_avg = np.mean(np.mean(syn_dat[:,x[0]:x[1],y[0]:y[1],0], axis = 1), axis = 1)
    
    #return syn_wav, syn_avg, obs_wav, obs_dat
    #wc_syn_idx = 26
    #wc_obs_idx = 10
    syn_ip = np.zeros((5,nw_obs, nx, ny))
    #dl = np.abs(obs_wav[wc_obs_idx]-syn_wav[wc_syn_idx])
    
    #w = [6,45]
    for s in range(4):
        syn_ip[s+1,:,:,:] = syn_dat[:,:,:,s]#np.interp(obs_wav, syn_wav[w[0]:w[1]]-dl, syn_dat[w[0]:w[1],i,j, s])
    if(calibrate):
        syn_avg_ip = np.mean(np.mean(syn_ip[1,:,:,:], axis = 1), axis = 1)
        
        syn_wav_cal, syn_dat_cal, factor, spec_fts, units = calib.spectrum(obs_wav #syn_wav[w[0]:w[1]]-dl
                                                                           , syn_ip[1,:,:,:] #fe_i[w[0]:w[1],:,:]
                                                                           , wave_idx=[0,50,123,-1]
                                                                           , extra_weight=40
                                                                           , atlas_range=0.5
                                                                           , spec_avg=syn_avg_ip #g[w[0]:w[1]]
                                                                           , calib_wave = True, cgs=True, verbose = True)
        syn_ip[1,:,:,:] = syn_dat_cal

    #stop()
    if(calibrate==1):
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_ip[0,:,i,j] = syn_wav_cal
        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_fe_6302_calibrated.fits'
    if(calibrate==0):
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_ip[0,:,i,j] = syn_wav
        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_fe_6302.fits'
        

# Ca 8542
if (line==8542):

    # definitions
    obs_file = mf.readfits(obs_dir+'b_obs8542_06_165307.fits')
    
    syn_dat = np.load(syn_dir + 'ca8_fullstokes_dat.npy')
    syn_wav = np.load(syn_dir + 'ca8_wav.npy')
    ca8_i = syn_dat[:,:,:,0]
    ca8_q = syn_dat[:,:,:,1]
    ca8_u = syn_dat[:,:,:,2]
    ca8_v = syn_dat[:,:,:,3]
    
    obs_wav = obs_file[0,:,0]
    obs_dat = obs_file[0:,:,:]
    obs_dat_avg = np.mean(obs_dat[1,:,:], axis = 1)
    nw_obs = len(obs_wav)
    
    nt = 1
    nx = syn_dat.shape[1]
    ny = syn_dat.shape[2]
    
    syn_avg = np.mean(np.mean(syn_dat[:,x[0]:x[1],y[0]:y[1],0], axis = 1), axis = 1)
    
    #return syn_wav, syn_avg, obs_wav, obs_dat
    wc_syn_idx = 26
    wc_obs_idx = 10
    syn_ip = np.zeros((5,len(syn_wav), nx, ny))
    #dl = np.abs(obs_wav[wc_obs_idx]-syn_wav[wc_syn_idx])

    if (calibrate==1):
        dl = np.abs(8542.09-syn_wav[wc_syn_idx])
    if(calibrate==0):
        dl = 0.
    
    w = [1,-1]
    for s in range(4):
        syn_ip[s+1,:,:,:] = syn_dat[:,:,:,s]
        
    if(calibrate):
        syn_avg_ip = np.interp(obs_wav, syn_wav[w[0]:w[1]]-dl, syn_avg[w[0]:w[1]])
        
        syn_wav_cal, syn_dat_cal, factor, spec_fts, units = calib.spectrum(syn_wav[w[0]:w[1]]-dl
                                                                           , ca8_i[w[0]:w[1],:,:]
                                                                           , wave_idx=[4,-5]
                                                                           , extra_weight=40
                                                                           , atlas_range=0.5
                                                                           , spec_avg=syn_avg[w[0]:w[1]] #g[w[0]:w[1]]
                                                                           #, calib_wave = True
                                                                           , cgs=True, verbose = True)
        syn_fin = np.zeros((5,len(syn_wav_cal), nx, ny))
        syn_fin[1,:,:,:] = syn_dat_cal
        syn_fin[2,:,:,:] = ca8_q[w[0]:w[1],:,:]
        syn_fin[3,:,:,:] = ca8_u[w[0]:w[1],:,:]
        syn_fin[4,:,:,:] = ca8_v[w[0]:w[1],:,:]
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_fin[0,:,i,j] = syn_wav_cal
        syn_ip = syn_fin
        #stop()

    if(calibrate==1):
        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_ca8_8542_calibrated.fits'
    if(calibrate==0):
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_ip[0,:,i,j] = syn_wav
        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_ca8_8542.fits'

# Ca II K
if (line==3950):
    
    # definitions
    obs_file = mf.readfits(obs_dir+'b_obs3950_06_165307.fits')
    
    syn_dat = np.load(syn_dir + 'cak_dat.npy')
    syn_wav = np.load(syn_dir + 'cak_wav.npy')
    cont_dat = np.load(syn_dir + 'cont_4000.11.npy')
    cont_wav = np.array([4000.11])#[4001.13951])
    cak_i = syn_dat[:,:,:]
    
    obs_wav = obs_file[0,:,0]
    obs_dat = obs_file[0:,:,:]
    obs_dat_avg = np.mean(obs_dat[1,:,:], axis = 1)
    nw_obs = len(obs_wav)
    
    nt = 1
    nx = syn_dat.shape[1]
    ny = syn_dat.shape[2]

    
    # Ca II K profile calibration
    # ===============
    syn_avg = np.mean(np.mean(syn_dat[:,x[0]:x[1],y[0]:y[1]], axis = 1), axis = 1)
    syn_avg_fov = np.mean(np.mean(syn_dat[:,:,:], axis = 1), axis = 1)

    #return syn_wav, syn_avg, obs_wav, obs_dat
    wc_syn_idx = 55
    wc_obs_idx = 10
    syn_ip = np.zeros((5,len(syn_wav), nx, ny))

    if(calibrate==1):
        dl = np.abs(3933.68-syn_wav[wc_syn_idx])#np.abs(obs_wav[wc_obs_idx]-syn_wav[wc_syn_idx])
    if(calibrate==0):
        dl = 0.
        
    w =[1,-1]# [np.argmin(np.abs(syn_wav-dl-obs_wav[0])),np.argmin(np.abs(syn_wav-dl-obs_wav[-2]))]
    
    syn_ip[1,:,:,:] = syn_dat[:,:,:]
        
    if(calibrate):
        syn_avg_ip = np.interp(obs_wav[:-1], syn_wav[w[0]:w[1]]-dl, syn_avg[w[0]:w[1]])
        
        syn_wav_cal, syn_dat_cal, factor, spec_fts, units = calib.spectrum(syn_wav[w[0]:w[1]]-dl #syn_wav[w[0]:w[1]]-dl
                                                                           , syn_dat[w[0]:w[1], :, :] #syn_dat[w[0]:w[1],:,:]
                                                                           , wave_idx=np.array([2,6,10,27,49,75,-7,-3])
                                                                           , extra_weight=40
                                                                           , atlas_range=0.5
                                                                           , spec_avg=syn_avg[w[0]:w[1]] #g[w[0]:w[1]]
                                                                           #, calib_wave = True
                                                                           , cgs=True, verbose = True)

        print('press c to proceed to calibrating the continuum at 4000.')
        stop()
        plt.close("all")
        
        # continuum Calibration
        cont_avg = np.array([np.mean(cont_dat[x[0]:x[1],y[0]:y[1]])])
        cont_wav_cal, cont_dat_cal, factor_cont, spec_fts_cont, units = calib.spectrum(cont_wav #syn_wav[w[0]:w[1]]-dl
                                                                                       , cont_dat #syn_dat[w[0]:w[1],:,:]
                                                                                       , wave_idx=np.array([0])
                                                                                       , extra_weight=40
                                                                                       , atlas_range=0.5
                                                                                       , spec_avg= cont_avg#g[w[0]:w[1]]
                                                                                       #, calib_wave = True
                                                                                       , cgs=True, verbose = True)
        
        #stop()
    
        syn_fin = np.zeros((5,len(syn_wav_cal)+1, nx, ny))
        syn_fin[1,-1,:,:] = cont_dat_cal
        syn_fin[1,:-1,:,:] = syn_dat_cal
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_fin[0,:-1,i,j] = syn_wav_cal
                syn_fin[0,-1,:,:] = cont_wav_cal
        syn_ip = syn_fin

    if(calibrate==1):
        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_cak_3950_calibrated.fits'    
    #stop()
    if(calibrate==0):
        for i in range(nx):                                                           
            for  j in range(ny):
                syn_ip[0,:,i,j] = syn_wav

        filename = '/home/seki2695/INV/stic/III/calibrated_cubes/MuRAMsyn_cak_3950.fits'

    
if (save_fits):
    mf.writefits(filename, syn_ip)
    print('.fits result is saved in '+filename)

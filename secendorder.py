from read_deal import *
from astropy.io import fits
from scipy import interpolate
from output_plot import t_iue
import numpy as np
import matplotlib.pyplot as plt
import os

def testFlux():
    juno = fits.open('/Users/zexixing/Research/swiftASTER/docs/juno/00091026003_2order.pha')
    data = juno[2].data
    rate = data['NETRATE']
    flux = data['FLUX']
    wave = data['LAMBDA']
    binwidth = data['BINWIDTH']
    arf = readarf(order=1)
    ea = arf(wave)
    ea0 = data['EFFAREA1']
    h_planck = 6.626e-27  # erg/s
    lightspeed = 2.9979e10  # cm/sec
    h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
    hnu = h_c_ang/(wave)
    flux_pre = hnu*rate/(ea*binwidth)
    fcoi = data['SP1_COIF']
    bgcoi = data['BG1_COIF']
    # plot
    plt.plot(wave,fcoi,'r-')
    plt.plot(wave,bgcoi,'b-')
    plt.show() #?????

def count2flux():
    pass

def remove2nd(aster,data_name,output=False):
    # readin data
    rela_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    data_path = rela_path + data_name
    data = fits.open(data_path)[2].data
    pix = data['PIXNO']
    wave = data['LAMBDA']
    netrate = data['NETRATE']
    bkgrate = data['BGRATE1']
    flux_uvotpy = data['FLUX']
    apercorr = data['APERCOR1']
    exposure = data['EXPOSURE']
    fcoi = data['SP1_COIF']
    bgcoi = data['BG1_COIF']
    senscorr = data['SENSCORR']
    binwidth = data['BINWIDTH']
    ea = data['EFFAREA1']
    wave2 = data['LAMBDA2']
    binwidth2 = data['BINWIDT2']
    ea2 = data['EFFAREA2']
    wave2[wave2==0]=np.NaN
    ea2_func = readarf(order=2)
    ea2_ = ea2_func(wave2)
    # calibrated count rates to newly calibrated count rates (src, bkg)
    # calibration: aperture; sensitivity; coincidence
    netrate = (netrate/fcoi + bkgrate/bgcoi)*fcoi - bkgrate
    # rate to flux
    h_planck = 6.626e-27  # erg/s
    lightspeed = 2.9979e10  # cm/sec
    h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
    hnu = h_c_ang/(wave)
    flux_recalib = hnu*netrate/(ea*binwidth)
    # remove 2nd order
    wave_flux = interpolate.interp1d(wave,flux_recalib,bounds_error=False,fill_value=np.NaN)
    hnu2 = h_c_ang/(wave2)
    netrate2 = wave_flux(wave2)*ea2_*binwidth2/hnu2
    netrate_corr2nd = netrate - netrate2
    flux_corr2nd = hnu*netrate_corr2nd/(ea*binwidth)
    plt.plot(pix,flux_recalib,label='1st order')
    plt.plot(pix,wave_flux(wave2),label='2nd order')
    plt.xlabel('pixel')
    plt.ylabel('erg s-1 cm-2 A-1')
    plt.ylim(-0.25e-12,1e-12)
    plt.legend()
    plt.show()
    # readin IUE
    data_iue = t_iue(['LWR01896LL.FITS','LWR05678LL.FITS','LWR05679LL.FITS','LWR05690LL.FITS','LWR06487LL.FITS'],plot=False)
    data_path_iue = rela_path + aster+'_raw_for_zexi.txt'
    #data_iue = readinSpec(data_path_iue, 'txt', 0, 1, sep='')
    wave_iue_ = np.array(data_iue['wave'])
    flux_iue = np.array(data_iue['flux'])
    wave_iue = wave_iue_[wave_iue_<3200]
    flux_iue = flux_iue[wave_iue_<3200]
    # remove 2nd order by IUE
    wave_flux_iue = interpolate.interp1d(wave_iue,flux_iue,bounds_error=False,fill_value=np.NaN)
    #----- wave limit
    netrate2_iue = wave_flux_iue(wave2)*ea2_*binwidth2/hnu2
    netrate_corr2nd_iue = netrate - netrate2_iue
    flux_corr2nd_iue = hnu*netrate_corr2nd_iue/(ea*binwidth)
    # plot
    fig = plt.figure()
    #plt.plot(wave,flux_uvotpy,label='uvotpy')
    plt.plot(wave,flux_recalib)#,label='recalibrated')
    plt.plot(wave,flux_corr2nd,label='remove 2nd order')
    plt.plot(wave_iue,flux_iue,label='IUE')
    #plt.plot(wave,flux_corr2nd_iue,label='remove 2nd order by IUE')
    plt.legend()
    plt.ylim(-0.25e-12,1e-12)
    plt.xlabel('wavelength')
    plt.ylabel('erg s-1 cm-2 A-1')
    plt.show()
    fig = plt.figure()
    plt.plot(wave[wave>2000],((flux_recalib-flux_corr2nd)/flux_recalib)[wave>2000]*100,label='remove 2nd order by obs')
    plt.plot(wave[wave>2000],((flux_recalib-flux_corr2nd_iue)/flux_recalib)[wave>2000]*100,label='remove 2nd order by iue')
    plt.ylim(-2.5,10)
    plt.ylabel('%')
    plt.xlabel('wavelength')
    plt.legend()
    plt.show()

#remove2nd('juno','00091026003_2order.pha')

mod_dict = {'wd0320':'wd0320_539_mod_001.fits',#'wd0320_539_stis_005.fits',
           'wd1057':'wd1057_719_mod_006.fits',
           'wd1657':'wd1657_343_mod_007.fits',
           'gd153':'gd153_mod_011.fits',
           'p177':'p177d_stisnic_008.fits',#'p177d_mod_003.fits',
           'p41':'p041c_stisnic_008.fits',#'p041c_001.fits'
           'agk81':'agk_81d266_stisnic_007.fits',
           'bd25':'bd_25d4655_002.fits',}

def teststar(star,obsid,norm=False):
    # readin mod
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
    name_mod = mod_dict[star]
    path_mod = docsdir + name_mod
    mod = fits.open(path_mod)
    data_mod = mod[1].data
    flux_mod = data_mod['FLUX']
    wave_mod = data_mod['WAVELENGTH']
    mod = interpolate.interp1d(wave_mod,flux_mod,bounds_error=False,fill_value=np.NaN)
    # readin obs
    if norm == False:
        path_obs = docsdir+obsid+'_1_default.pha'
        path_bkg = docsdir+obsid+'_1_default_back.pha'
    else:
        path_obs = docsdir+obsid+'_1_default_norm.pha'
        path_bkg = docsdir+obsid+'_1_default_norm_back.pha' 
    obs = fits.open(path_obs)
    bkg = fits.open(path_bkg)
    data_obs = obs[2].data
    counts = obs[1].data['COUNTS'][::-1]
    bkgcounts = bkg[1].data['BKG_COUNTS'][::-1]
    pixel = data_obs['PIXNO']
    wave = data_obs['LAMBDA']
    flux = data_obs['FLUX']
    netrate = data_obs['NETRATE']
    bgrate = data_obs['BGRATE1']
    aper = data_obs['APERCOR1']
    coif = data_obs['SP1_COIF']
    coibkg = data_obs['BG1_COIF']
    sens = data_obs['SENSCORR']
    exp = data_obs['EXPOSURE']
    eff = data_obs['EFFAREA1']
    binwidth = data_obs['BINWIDTH']
    wave2 = data_obs['LAMBDA2']
    binwidth2 = data_obs['BINWIDT2']
    #wave2[wave2==0]=np.NaN
    pixel2 = data_obs['PIXNO2']
    #pixel2[pixel2==0]=np.NaN
    sens2 = data_obs['SENSCOR2']
    coif2 = data_obs['SP2_COI']
    aper2 = data_obs['APERCOR2']
    exp2 = data_obs['EXPOSURE2']
    #-------
    h_planck = 6.626e-27  # erg/s
    lightspeed = 2.9979e10  # cm/sec
    h_c_ang = h_planck * lightspeed * 1e8 # ergs-Angstrom
    hnu = h_c_ang/(wave) 
    #--------
    if norm == False:
        swifttime = obs[1].header['TSTART']
        fscale =(swifttime-126230400.000) / (12.6*365.26*86400)
    else:
        fscale = 1
    # pred
    flux_mod = mod(wave)
    counts_pred = flux_mod*binwidth*eff*exp/(hnu*sens*fscale*aper) #/coif #+ bkgcounts
    wave_pred = wave
    #pixel_pred = pixel
    #plot
    #plt.plot(wave,counts,'b-')
    #plt.plot(wave,counts_pred,'r-')
    #plt.show()
    # 2nd order
    diff = pixel2[0]-pixel[0]
    eff2_func = readarf(order=2)#,smooth=False)
    eff2 = eff2_func(wave2)
    flux2_mod = mod(wave2)
    counts2_pred = flux2_mod*binwidth2*eff2*exp2/(hnu*sens2*fscale*aper2)  #/coif2
    counts2_pred = np.append(np.zeros(diff),counts2_pred[:-diff])
    #pixel2_pred = np.append(np.arange(pixel[0],pixel2[0]),pixel2[:-diff])
    #wave2_pred = np.append(np.zeros(diff),wave2[:-diff])
    counts_pred = np.nan_to_num(counts_pred)
    counts2_pred = np.nan_to_num(counts2_pred)
    counts_total_pred = (counts2_pred+counts_pred)/coif

    eff22_func = readarf(order=2,smooth=False)
    eff22 = eff22_func(wave2)
    counts22_pred = flux2_mod*binwidth2*eff22*exp2/(hnu*sens2*fscale*aper2)*0.8 #/coif2
    counts22_pred = np.append(np.zeros(diff),counts22_pred[:-diff])
    plt.plot(wave,(counts22_pred+counts_pred)/coif,'r-',alpha=0.5)
    plt.plot(wave,counts22_pred,'k-',alpha=0.2)

    plt.plot(wave,counts_total_pred,'r-')
    plt.plot(wave,counts-bkgcounts,'b-')
    plt.plot(wave,counts_pred,'k-',alpha=0.2)
    plt.plot(wave,counts2_pred,'k-',alpha=0.2)
    #plt.plot(wave,bkgcounts,'k-',alpha=0.2)
    plt.xlabel('wavelength')
    plt.ylabel('counts')
    plt.title(star+' '+obsid)
    plt.ylim(0)
    plt.show()

    flux_meas = (counts-bkgcounts)*(hnu*sens*fscale*aper*coif)/(binwidth*eff*exp)
    plt.plot(wave,flux_mod,'r-')
    plt.plot(wave,flux_meas,'b-')
    plt.ylim(-0.2e-13,1e-13)
    plt.show()
    
    refl = np.nan_to_num(flux_meas/flux_mod,nan=-100)
    wave = wave[refl!=-100]
    refl = refl[refl!=-100]
    spec_dict = binSpec(wave, refl, binGenerate(wave, 20), sigma=None)
    plt.plot(spec_dict['wave'],spec_dict['flux'])
    plt.ylim(0.25,1.75)
    plt.hlines(1,xmin=min(wave),xmax=max(wave))
    plt.show()

    #plt.plot(wave,eff,'r-')
    #plt.plot(wave2,eff2,'b-')
    #plt.plot(wave2,ea2_func(wave2),'b-',alpha=0.3)
    #plt.show()

def analogred():
    #mod = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    wave_norm = 3800
    mod=fits.open('/Users/zexixing/Research/swiftASTER/docs/standard/bd124536/gs_126.fits')
    mod = {'wave':mod[1].data['WAVELENGTH'],'flux':mod[1].data['FLUX']}
    mod = normRefl(mod, wave_norm=wave_norm)
    analog = getSolarFlux('00091706004.pha', if_smooth=True)
    analog = normRefl(analog, wave_norm=wave_norm)
    #mod_f = interpolate.interp1d(mod['wave'],mod['flux'])
    #mod['wave'] = analog['wave']
    #mod['flux'] = mod_f(mod['wave'])
    #refl = analog['flux']/mod['flux']
    #refl = normRefl({'wave':mod['wave'],'flux':refl}, wave_norm=3800)
    #plt.plot(refl['wave'],refl['flux'])
    #plt.ylim(0.25,1.75)
    plt.plot(mod['wave'],mod['flux'])
    plt.plot(analog['wave'],analog['flux'])
    plt.show()


analogred()
    
#teststar('wd0320','00054250027', norm=True)
#teststar('wd0320','00054250023')
#teststar('agk81','00057530014', norm=True)
#teststar('agk81','00057530004')
#teststar('p177','00056760024')
#teststar('p41','00057967002')
#teststar('p41','00057955002')
#teststar('wd1057','00055205001')
#teststar('wd1057','00055211004')
#teststar('wd1057','00055216001')
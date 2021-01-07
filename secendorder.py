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



remove2nd('juno','00091026003_2order.pha')
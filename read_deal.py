# plot spectra
import pandas as pd
import numpy as np
from tools import *
from scipy import interpolate
from astropy.io import fits
import matplotlib.pyplot as plt

remove_list = [#'00091503001', #flora -- star
                '00091223001','00091225001','00091229001', #pallas -- star(compact)
                '00091198002', #'00091022002','00091197002', vesta -- star
                '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                '00091540001','00091538001', #nysa -- star(compact)
                '00091591001', #themis -- wrong position
                '00091207002', #juno -- star in slit
                '00091237001', #dembowska -- star in slit
                '00091268001', '00091501001', '00091503001', '00091507001', #flora -- star in slit
                '00091559001', #hygiea -- star in slit
                '00091519001', '00091521001', '00091523001', '00091525001', #iris -- star in slit
                '00091598001', '00091593001', '00091595001', #themis -- star in slit
                ]

def readinSpec(spec_path, filetype, wave_index, flux_index,
               sigma_index=None, flag_index=None,
               ext=1, sep=' ', skiprows=0):
    # initialize spec_dict
    spec_dict = {}
    # from fits
    if filetype == 'fits':
        hdul = fits.open(spec_path)
        hdu = hdul[ext].data
        #if wave_index in hdu.names and flux_index in hdu.names:
        spec_dict['wave'] = hdu[wave_index]
        spec_dict['flux'] = hdu[flux_index]
        #else:
        #    raise ValueError('wrong wave/flux keywords')
        if sigma_index:
            spec_dict['sigma'] = hdu[sigma_index]
        if flag_index:
            spec_dict['flag'] = hdu[flag_index]
    # from csv (no words in the table)
    elif filetype == 'csv':
        data = np.loadtxt(spec_path, delimiter=sep, skiprows=skiprows)
        if isinstance(wave_index,int) and isinstance(flux_index ,int):
            spec_dict['wave'] = data[:,wave_index]
            spec_dict['flux'] = data[:,flux_index]
        else:
            raise TypeError('wrong wave/flux keywords')
        if isinstance(sigma_index,int):
            spec_dict['sigma'] = data[:,sigma_index]
        if isinstance(flag_index,int):
            spec_dict['flag'] = data[:,flag_index]
    # from txt (some words in the table)
    elif filetype == 'txt':
        f = open(spec_path,'rb')
        wave, flux, sigma, flag = [], [], [], []
        for line in f.readlines()[skiprows:]:
            line = line.strip()
            #line = line.strip('\n')
            line = line.split()
            if isinstance(wave_index,int) and isinstance(flux_index ,int):
                wave.append(float(line[wave_index]))
                flux.append(float(line[flux_index]))
            else:
                raise TypeError('wrong wave/flux keywords')
            if isinstance(sigma_index,int):
                sigma.append(float(line[sigma_index]))
            if isinstance(flag_index,int):
                flag.append(float(line[flag_index]))
        spec_dict['wave'] = wave
        spec_dict['flux'] = flux
        if isinstance(sigma_index,int):
            spec_dict['sigma'] = sigma
        if isinstance(flag_index,int):
            spec_dict['flag'] = flag
    else:
        pass
    return spec_dict

def binGenerate(wave, bin_step):
    bins = np.arange(wave[0],wave[-1]+bin_step,bin_step)
    return bins

def binSpec(wave, flux, bins, sigma=None, flag=None):
    spec_dict = {}
    wave = np.array(wave)
    flux = np.array(flux)
    bins = np.array(bins)
    if sigma is not None:
        sigma = np.array(sigma)
    if flag is not None:
        flag = np.array(flag)
        wave = wave[flag==0]
        flux = flux[flag==0]
        if sigma is not None:
            sigma = sigma[flag==0]
    groupIdx = np.digitize(wave, bins)
    wave_means = [wave[groupIdx == i].mean() for i in range(1, len(bins))]
    flux_means = [flux[groupIdx == i].mean() for i in range(1, len(bins))]
    spec_dict['wave'] = np.array([i for i in wave_means if not np.isnan(i)])
    spec_dict['flux'] = np.array([i for i in flux_means if not np.isnan(i)])
    if sigma is not None:
        sigma_means = [np.power(sigma[groupIdx == i],2) for i in range(1, len(bins))]
        sigma_means = [np.sqrt(sum(i))/len(i) for i in sigma_means]
        spec_dict['sigma'] = np.array([i for i in sigma_means if not np.isnan(i)])
    return spec_dict

def mergeSpec(wave_list, flux_list, sigma_list=None, flag_list=None):
    spec_dict = {}
    wave = [item for sublist in wave_list for item in sublist]
    flux = [item for sublist in flux_list for item in sublist]
    spec_dict['wave'] = wave
    spec_dict['flux'] = flux
    if sigma_list is not None:
        sigma = [item for sublist in sigma_list for item in sublist]
        spec_dict['sigma'] = sigma
    if flag_list is not None:
        flag = [item for sublist in flag_list for item in sublist]
        spec_dict['flag'] = flag
    return spec_dict

def normRefl(spec_dict, wave_norm=3800):
    f = interpolate.interp1d(spec_dict['wave'],spec_dict['flux'])
    norm_f = f(wave_norm)
    spec_dict['flux'] = spec_dict['flux']/norm_f
    if 'sigma' in spec_dict.keys():
        #s = interpolate.interp1d(spec_dict['wave'],spec_dict['sigma'])
        #norm_s = s(wave_norm)
        spec_dict['sigma'] = spec_dict['sigma']/norm_f
    return spec_dict

def unitConv(spec_dict, unit={'energy':'erg/s','area':'cm-2'}):
    # spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], 10), sigma=spec_dict['sigma'], flag=spec_dict['flag'])
    if unit['energy'] == 'W':
        spec_dict['flux'] = np.array(spec_dict['flux'])*1e-7
        spec_dict['sigma'] = np.array(spec_dict['sigma'])*1e-7
    if unit['area'] == 'm-2':
        spec_dict['flux'] = np.array(spec_dict['flux'])*1e4
        spec_dict['sigma'] = np.array(spec_dict['sigma'])*1e4
    return spec_dict

def interpSpec(spec_dict, wave_bins=binGenerate((1660,6800), 20)):
    # TODO: sigma & flux 
    f = interpolate.interp1d(spec_dict['wave'],spec_dict['flux'],fill_value="extrapolate")
    spec_dict['wave'] = wave_bins
    spec_dict['flux'] = f(wave_bins)
    return spec_dict

def readarf(order=1, smooth=True, iferr=False, filt='uvgrism'):
    if filt == 'uvgrism':
        if order == 1:
            if smooth==True:
                arf_name = 'swugu0160_20041120v105.arf'
                extname = 'SPECRESPUGRISM0160_AX1040AY1030_O1'
            else:
                arf_name = 'swugu0160_1_20041120v999.arf'
                extname = 'SPECRESPUGRISM160'
        elif order == 2:
            extname = 'SPECRESP0160GRISM2NDORDER'
            if smooth==True:
                arf_name = 'swugu0160_ax1080ay1000_dx100dy100_o2_20041120v001.arf'
            else:
                arf_name = 'swugu0160_2_20041120v999.arf'
        rela_path = '/Users/zexixing/pymodules/uvotpy-2.3.1/uvotpy/calfiles/'
    else:
        rela_path = '/Users/zexixing/Research/swiftASTER/docs/auxil/'
        arf_name = 'arf_'+filt+'.fits'
        extname=1
    path = rela_path+arf_name
    arf = fits.open(path)[extname].data
    wave = (arf['WAVE_MIN']+arf['WAVE_MAX'])/2
    ea = arf['SPECRESP']
    f_arf = interpolate.interp1d(wave,ea,bounds_error=False,fill_value=np.NaN)
    if iferr == False:
        return f_arf
    else:
        eaerr = arf['SPRS_ERR']
        f_err = interpolate.interp1d(wave,eaerr,bounds_error=False,fill_value=np.NaN)
        return f_arf, f_err


def remove2ndFlux(aster,data_name):
    # readin data
    rela_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    data_path = rela_path + data_name
    data = fits.open(data_path)[2].data
    wave = data['LAMBDA']
    netrate = data['NETRATE']
    bkgrate = data['BGRATE1']
    flux_uvotpy = data['FLUX']
    flux_err = data['FLUXERR']
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
    flux_uvotpy = flux_uvotpy[wave!=0]
    flux_recalib = flux_recalib[wave!=0]
    flux_corr2nd = flux_corr2nd[wave!=0]
    wave = wave[wave!=0]
    dict = {'wave':wave,
            'flux_uvotpy':flux_uvotpy,
            'sigma':flux_err,
            'flux_recalib':flux_recalib,
            'flux_corr2nd':flux_corr2nd}
    return dict

def coadd_(spec_list,aster='.',binbox=10):
    wave_list = []
    flux_list = []
    sigma_list = []
    flag_list = []
    if type(spec_list[0]) == str:
        for name in spec_list:
            spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+name, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
            wave_list.append(spec_dict['wave'])
            flux_list.append(spec_dict['flux'])
            sigma_list.append(spec_dict['sigma'])
            flag_list.append(spec_dict['flag'])
        spec_dict = mergeSpec(wave_list, flux_list, sigma_list=sigma_list, flag_list=flag_list)
        spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], binbox), 
                            sigma=spec_dict['sigma'], flag=spec_dict['flag'])
    else:
        for spec_dict in spec_list:
            wave_list.append(spec_dict['wave'])
            flux_list.append(spec_dict['flux'])
            if 'sigma' in spec_dict.keys():
                sigma_list.append(spec_dict['sigma'])
            else:
                sigma_list = None
            #flag_list.append(spec_dict['flag'])
        spec_dict = mergeSpec(wave_list, flux_list, sigma_list=sigma_list, flag_list=None)
        if 'sigma' in spec_dict.keys():
            sigma = spec_dict['sigma']
        else:
            sigma = None
        spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], binbox), 
                            sigma=sigma, flag=None)
    return spec_dict

def coadd(spec_list,aster='.',remove2nd=False):
    wave_list = []
    flux_list = []
    sigma_list = []
    flag_list = []
    if type(spec_list[0]) == str:
        for name in spec_list:
            spec_dict = remove2ndFlux(aster,name)
            wave_list.append(spec_dict['wave'])
            if remove2nd:
                flux_list.append(spec_dict['flux_corr2nd'])
            else:
                flux_list.append(spec_dict['flux_uvotpy'])
        spec_dict = mergeSpec(wave_list, flux_list)
        spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], 10))
    else:
        for spec_dict in spec_list:
            wave_list.append(spec_dict['wave'])
            if remove2nd:
                flux_list.append(spec_dict['flux_corr2nd'])
            else:
                flux_list.append(spec_dict['flux_uvotpy'])
            #sigma_list.append(spec_dict['sigma'])
            #flag_list.append(spec_dict['flag'])
        spec_dict = mergeSpec(wave_list, flux_list, sigma_list=None, flag_list=None)
        spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], 10), 
                            sigma=None, flag=None)
    return spec_dict

def lidird_solar(date, file_path):
    data = pd.read_csv(file_path, skiprows=1, names=['day','wave','flux'])
    #data['day']=data['day'].astype(int)
    data = dict(list(data.groupby('day')))
    wave = np.array(data[date]['wave'])*10
    flux = np.array(data[date]['flux'])/10
    f = interpolate.interp1d(wave,flux,fill_value="extrapolate")
    integ = sum(flux)*(wave[2]-wave[1])
    return f, integ
    #% example
    #file_path = '/Users/zexixing/Research/swiftASTER/docs/omi_ssi.csv'
    #date=3919
    #f1 = lidird_solar(date, file_path)[0]
    #f2 = lidird_solar(date+366+366, file_path)[0]
    #f3 = lidird_solar(date+365*3+366, file_path)[0]
    #f4 = lidird_solar(date+365*4+366*2, file_path)[0]
    #f5 = lidird_solar(date+365*6+366*2, file_path)[0]
    #wave = np.arange(2650,5000+10,10)
    #plt.plot(wave,f1(wave),label='2006-7-2')
    #plt.plot(wave,f2(wave),label='2008-7-2')
    #plt.plot(wave,f3(wave),label='2010-7-2')
    #plt.plot(wave,f4(wave),label='2012-7-2')
    #plt.plot(wave,f5(wave),label='2014-7-2')
    #plt.legend()
    #plt.show()
    
def ugrismlist(aster,remove=False, filt='uvgrism'):
    datadir = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'#+obsid+'/uvot/image'
    namelist_ = os.listdir(datadir)
    namelist_ = [i for i in namelist_ if i[:3]=='000']
    if filt=='uvgrism':
        namelist = []
        for obsid in namelist_:
            obsdir = datadir+obsid+'/uvot/image'
            filelist = os.listdir(obsdir)
            if ('sw'+obsid+'ugu_sk.img.gz' in filelist) or \
                ('sw'+obsid+'ugu_sk.img' in filelist):
                namelist.append(obsid)
    elif filt=='allfilter':
        namelist = [i[:11] for i in namelist_]
    else:
        namelist = []
        for obsid in namelist_:
            obsdir = datadir+obsid+'/uvot/image'
            filelist = os.listdir(obsdir)
            if 'sw'+obsid+filt+'_rw.img.gz' in filelist:
                namelist.append(obsid)
    if remove != False: 
        #remove_list = ['00091503001', #flora -- star
        #               '00091223001','00091225001','00091229001', #pallas -- star(compact)
        #               '00091198002', #'00091022002','00091197002', vesta -- star
        #               '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
        #               '00091540001','00091538001', #nysa -- star(compact)
        #               '00091591001', #themis -- wrong position
        #               ]
        for obsid_rm in remove:
            try:    namelist.remove(obsid_rm)
            except:    pass
    return namelist

def colorFader(mix=0,c1='blue',c2='red'): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    import matplotlib as mpl
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def getSolarFlux(solar_spec_name, if_smooth=True):
    # read solar spec
    if solar_spec_name[-4:] == '.pha':
        solar_path = '/Users/zexixing/Research/swiftASTER/docs/'+solar_spec_name
        solar_spec_dict = readinSpec(solar_path, 
                                    'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                    flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
        if if_smooth == True:
            solar_spec_dict=coadd_([solar_spec_dict])
        else:
            pass
    elif solar_spec_name == 'sun_ref_colina96.asc.txt':
        solar_spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/sun_ref_colina96.asc.txt', 'csv', 0, 1,
                                     sigma_index=None, flag_index=None,
                                     ext=1, sep='    ', skiprows=0)
    else:
        solar_spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/'+solar_spec_name, 'csv', 0, 1,
                                     sigma_index=None, flag_index=None,
                                     ext=1, sep=',', skiprows=0)     
    #solar_spec_dict = unitConv(solar_spec_dict, unit={'energy':'erg/s','area':'cm-2'})
    return solar_spec_dict

def getRefl(aster_spec_name, solar_spec_dict, aster='.'):
    from tools import error_prop
    # TODO: sigma & flux
    # read the asteroid's spec
    if type(aster_spec_name) == str:
        aster_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+aster_spec_name
        aster_spec_dict = readinSpec(aster_path, 
                                    'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                    flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    else:
        aster_spec_dict = aster_spec_name
    #aster_spec_dict = unitConv(aster_spec_dict, unit={'energy':'erg/s','area':'cm-2'})
    f = interpolate.interp1d(solar_spec_dict['wave'],solar_spec_dict['flux'],fill_value="extrapolate")
    refl = aster_spec_dict['flux']/f(aster_spec_dict['wave'])
    wave = aster_spec_dict['wave']
    #aster_spec_dict = interpSpec(aster_spec_dict)
    #solar_spec_dict = interpSpec(solar_spec_dict)
    #refl = aster_spec_dict['flux']/solar_spec_dict['flux']
    #wave = aster_spec_dict['wave']
    # -----------error
    refl_dict = {'wave':wave,'flux':refl}
    if 'sigma' in solar_spec_dict.keys():
        e = interpolate.interp1d(solar_spec_dict['wave'],solar_spec_dict['sigma'],fill_value="extrapolate")
        sigma = error_prop('div', aster_spec_dict['flux'], aster_spec_dict['sigma'], f(aster_spec_dict['wave']), e(aster_spec_dict['wave']))
        refl_dict['sigma'] = sigma
    else:
        sigma = aster_spec_dict['sigma']/f(aster_spec_dict['wave'])
        refl_dict['sigma'] = sigma
    refl_dict = normRefl(refl_dict)
    if type(aster_spec_name) == str:
        return refl_dict
    else:
        return refl_dict
        #return refl_dict['wave'], refl_dict['flux']

def redden(wave, flux, reddening, l0=5500, l1=3000):
    x = (1-reddening/(20000/(l0-l1))*1)/(1+reddening/(20000/(l0-l1))*1)
    f = interpolate.interp1d((3000,5500),(x,1),fill_value="extrapolate")
    flux = f(wave)*flux
    return wave, flux

def mean_aster_flux(aster,adjname = '_smeargauss',binbox=10):
    idlist = ugrismlist(aster)
    #remove_list = [#'00091503001', #flora -- star
    #               '00091223001','00091225001','00091229001', #pallas -- star(compact)
    #               '00091198002', #'00091022002','00091197002', vesta -- star
    #               '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
    #               '00091540001','00091538001', #nysa -- star(compact)
    #               '00091591001', #themis -- wrong position
    #               '00091207002', #juno -- star in slit
    #               '00091237001', #dembowska -- star in slit
    #               '00091268001', '00091501001', '00091503001', '00091507001', #flora -- star in slit
    #               '00091559001', #hygiea -- star in slit
    #               '00091519001', '00091521001', '00091523001', '00091525001', #iris -- star in slit
    #               '00091598001', '00091593001', '00091595001', #themis -- star in slit
    #               ]
    for obsid_rm in remove_list:
        try:    idlist.remove(obsid_rm)
        except:    pass
    # loop
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    refl_list = []
    flux_list = []
    for obsid in idlist:
        # flux
        aster_spec_name = obsid+adjname+'.pha'
        spec_path= docsdir+aster_spec_name
        hdul = fits.open(spec_path)
        hdu = hdul[2].data
        flux_dict = {'wave':hdu['LAMBDA'],'flux':hdu['FLUX'],'sigma':hdu['FLUXERR']}
        flux_list.append(flux_dict)
    flux_dict = coadd_(flux_list,aster,binbox=binbox)
    return flux_dict

def solar_flux(sun, binbox=20):
    if sun == 'analog':
        solar_spec_name = '00091706004.pha'
    elif sun == 'colina96':
        solar_spec_name = 'sun_ref_colina96.asc.txt'
    else:
        raise ValueError('Wrong name of the solar flux spectrum!')
    solar_spec_dict = getSolarFlux(solar_spec_name, if_smooth=False)
    if binbox == False:
        return solar_spec_dict
    else:
        solar_spec_dict = coadd_([solar_spec_dict],binbox=binbox)
        return solar_spec_dict

def read_fits(aster, obsid, adjname):
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    aster_spec_name = obsid+adjname+'.pha'
    spec_path= docsdir+aster_spec_name
    hdul = fits.open(spec_path)
    return hdul

def aster_flux(aster, obsid, adjname, binbox=20):
    hdul = read_fits(aster, obsid, adjname)
    hdu = hdul[2].data
    flux_dict = {'wave':hdu['LAMBDA'],'flux':hdu['FLUX'],'sigma':hdu['FLUXERR']}
    if binbox == False:
        return flux_dict
    else:
        flux_dict = coadd_([flux_dict],aster,binbox=binbox)
        return flux_dict
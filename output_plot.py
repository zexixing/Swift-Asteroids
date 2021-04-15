from read_deal import *
import os
import pandas as pd

def plotSepc(wave, flux, color='b', linestyle='-',alpha=1, linewidth=0.7, label=None,
             sigma=None, flag=None):
    wave = np.array(wave)
    flux = np.array(flux)
    if sigma is not None:
        sigma = np.array(sigma)
    if flag is not None:
        flag = np.array(flag)
        wave = wave[flag==0]
        flux = flux[flag==0]
        if sigma is not None:
            sigma = sigma[flag==0]
    plt.plot(wave, flux, color=color, linestyle=linestyle, alpha=alpha,label=label,linewidth=linewidth)
    if label:
        plt.legend()
        pass
    if sigma is not None:
        pass #~TODO:
        #plt.fill_between(wave, flux-sigma, flux+sigma, color=color, alpha=0.6)
    # units

def printHdr(spec_path, ext, temp_path='/Users/zexixing/Downloads/hdr_temp.txt'):
    hdul = fits.open(spec_path)
    hdr = hdul[ext].header
    print(hdr)

def t_uvot(name,aster,color,linestyle='-',alpha=1,legend=None,unit={'energy':'erg/s','area':'cm-2'},remove2nd=False,show=1):
    #path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+name
    #spec_dict = readinSpec(path, 
    #                    'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
    #                    flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    #wave=spec_dict['wave']
    #flux=spec_dict['flux']
    #sigma=spec_dict['sigma']
    #flag=spec_dict['flag']
    spec_dict = remove2ndFlux(aster,name)
    wave=spec_dict['wave']
    if remove2nd == False:
        flux = spec_dict['flux_uvotpy']
    else:
        flux = spec_dict['flux_corr2nd']
    err = spec_dict['sigma']
    #spec_dict = binSpec(wave, flux, binGenerate(wave, 10), sigma=sigma, flag=flag)
    #wave = np.array(wave)
    #flux = np.array(flux)
    f = interpolate.interp1d(wave,flux)
    #err = np.array(spec_dict['sigma'])
    if unit['energy'] == 'W':
        flux = flux*1e-7
        err = err*1e-7
    if unit['area'] == 'm-2':
        flux = flux*1e4
        err = err*1e4
    plotSepc(wave, flux*1e12, 
             linestyle=linestyle,alpha=alpha,color=color, label=legend, sigma=err*1e12, flag=None)
    ylabel = 'flux ('+unit['energy']+' '+unit['area']+' A-1)'
    plt.ylabel(ylabel)

def t_iue(spec_list,plot=True):
    wave_list = []
    flux_list = []
    sigma_list = []
    flag_list = []
    for name in spec_list:
        spec_dict = readinSpec('/Users/zexixing/Dropbox/PDART Small Bodies UV Archive/IUE/data/raw/3 Juno/'+name, 
                              'fits', 'WAVELENGTH', 'FLUX', sigma_index='SIGMA', 
                              flag_index='QUALITY', ext=1, sep=' ', skiprows=0)
        wave_list.append(spec_dict['wave'])
        flux_list.append(spec_dict['flux'])
        sigma_list.append(spec_dict['sigma'])
        flag_list.append(spec_dict['flag'])
    spec_dict = mergeSpec(wave_list, flux_list, sigma_list=sigma_list, flag_list=flag_list)
    spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], 10), sigma=spec_dict['sigma'], flag=spec_dict['flag'])
    if plot == True:
        plotSepc(spec_dict['wave'], spec_dict['flux'], color='k', label='IUE',linestyle='-',sigma=spec_dict['sigma'], flag=None)
    else:
        return spec_dict

def writeFile(data_name,back_name, output_name):
    #----------write raw files----------
    # read in .pha and back.pha
    folder_path = '/Users/zexixing/Research/swiftASTER/docs/'
    data_hdul = fits.open(folder_path+data_name)
    data_dt = data_hdul[1].data
    ind = data_dt['CHANNEL']
    src = data_dt['COUNTS']
    src_err = data_dt['STAT_ERR']
    back_hdul = fits.open(folder_path+back_name)
    back_dt = back_hdul[1].data
    bkg = back_dt['BKG_COUNTS']
    bkg_err = back_dt['BKG_STAT_ERR']
    net = src - bkg
    net_err = error_prop('sub', src, src_err, bkg, bkg_err)
    # write files
    raw = np.array(list(zip(ind,src,src_err,bkg,bkg_err,net,net_err)))
    raw[:,0] = raw[:,0].astype(int)
    header = [np.array(['Pixel','Source','Delta Source','Background','Delta Background','Net','Delta Net']),
              np.array(['#','Counts','Counts','Counts','Counts','Counts','Counts'])]
    raw = pd.DataFrame(raw, columns=header)
    raw.Pixel = raw.Pixel.astype(int)
    raw.Source = raw.Source.astype(int)
    raw.Background = raw.Background.astype(int)
    raw.Net = raw.Net.astype(int)
    #raw.to_csv(folder_path+output_name,index=False)
    plt.show()

#writeFile('00091206002_1.pha','sw00091206002_1_back.pha', '00091206002.csv')

def compareAnalog():
    analog_list = [#'sw00091701002ugu_1ord_1_g_a.pha',
                   #'sw00091701002ugu_1ord_1_g_b.pha',
                   #'sw00091702002ugu_1ord_1_f.pha',
                   #'sw00091703002ugu_1ord_1_f.pha',
                   #'sw00091705002ugu_1ord_1_f.pha',
                   #'sw00091706004ugu_1ord_1_f.pha',
                   '00091706004_2order.pha']
                   #'sw00091738001ugu_1ord_1_f_a.pha',]
                   #'sw00091738001ugu_1ord_1_f_b.pha']
    wave_norm=3600
    color_list = ['gold','crimson', 'blueviolet', 'dodgerblue','aquamarine','greenyellow']#,'yellowgreen']
    linestyle_list = ['--','-.','--','-.','--',':']#,':']
    num = len(analog_list)
    for i in range(0,num):
        solar_path = '/Users/zexixing/Research/swiftASTER/docs/solar_analog/'+analog_list[i]
        spec_dict = readinSpec(solar_path, 
                               'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                               flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
        spec_dict = coadd_([spec_dict],binbox=20)
        spec_dict = normRefl(spec_dict,wave_norm=wave_norm)
        plotSepc(spec_dict['wave'], spec_dict['flux'], color='r', linestyle='-',alpha=1, linewidth=0.5,label='solar analog (Swift)',#analog_list[i][2:13],#color=color_list[i]
                 sigma=None, flag=None)
    # colina
    colina_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/sun_ref_colina96.asc.txt', 'csv', 0, 1,
                           sigma_index=None, flag_index=None,
                           ext=1, sep='    ', skiprows=0)
    colina_dict = coadd_([colina_dict],binbox=20)
    colina_dict = normRefl(colina_dict,wave_norm=wave_norm)
    colina_dict['wave'] = colina_dict['wave']
    plotSepc(colina_dict['wave'], colina_dict['flux'], color='b', linestyle='-',alpha=1, linewidth=0.5,label='solar model (colina96)',
             sigma=None, flag=None)
    #susim
    '''
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/susim_hires_norm.fits', 'fits', 0, 1,
                           sigma_index=None, flag_index=None,
                           ext=0, sep='    ', skiprows=0)
    spec_dict = normRefl({'wave':spec_dict['wave']*10,'flux':spec_dict['flux']},wave_norm=wave_norm)
    spec_dict = binSpec(spec_dict['wave'], np.array(spec_dict['flux']), binGenerate(spec_dict['wave'], 10), sigma=None)
    plotSepc(spec_dict['wave'], spec_dict['flux'], color='k', linestyle='-.',alpha=1, label='SUSIM',
             sigma=None, flag=None)
    
    # e490
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/sun_e490.csv', 'csv', 0, 1,
                           sigma_index=None, flag_index=None,
                           ext=1, sep=',', skiprows=1)
    spec_dict = normRefl({'wave':spec_dict['wave']*10000,'flux':spec_dict['flux']/10000},wave_norm=wave_norm)
    plotSepc(spec_dict['wave'], spec_dict['flux'], color='k', linestyle='--',alpha=1, label='E490',
             sigma=None, flag=None)
    # wehrli
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/sun_wehrli.csv', 'csv', 0, 1,
                           sigma_index=None, flag_index=None,
                           ext=1, sep=',', skiprows=1)
    spec_dict = normRefl({'wave':spec_dict['wave']*10,'flux':spec_dict['flux']/10},wave_norm=wave_norm)
    plotSepc(spec_dict['wave'], spec_dict['flux'], color='k', linestyle=':',alpha=1, label='wehrli',
             sigma=None, flag=None)
    '''
    aster='juno'
    idlist = ugrismlist(aster)
    remove_list = ['00091503001', #flora -- star
                   '00091223001','00091225001','00091229001', #pallas -- star(compact)
                   '00091198002', #'00091022002','00091197002', vesta -- star
                   '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                   '00091540001','00091538001', #nysa -- star(compact)
                   '00091591001', #themis -- wrong position
                   '00091207002', #juno -- star in slit
                   '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                   ]
    for obsid_rm in remove_list:
        try:    idlist.remove(obsid_rm)
        except:    pass
    # loop
    adjname = '_smeargauss'
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
    flux_dict = coadd_(flux_list,aster)
    flux_dict = normRefl(flux_dict, wave_norm=3600)
    plotSepc(flux_dict['wave'], flux_dict['flux']*1.45, color='k', label='Juno (Swift)',linestyle='-',sigma=None, flag=None)
    # refl - colina
    #refl_dict = getRefl(flux_dict, colina_dict)
    #refl_dict = normRefl(refl_dict, wave_norm=3800)
    #plotSepc(refl_dict['wave'], refl_dict['flux'], color='g', label='Refl',linestyle='-',sigma=None, flag=None)
    #f = interpolate.interp1d(refl_dict['wave'],refl_dict['flux'],fill_value="extrapolate")
    # refl - solar analog
    #refl_dict = getRefl(flux_dict, spec_dict)
    #refl_dict = normRefl(refl_dict, wave_norm=3800)
    #plotSepc(refl_dict['wave'], refl_dict['flux'], color='b', label='Refl',linestyle='-',sigma=None, flag=None)
    #f = interpolate.interp1d(refl_dict['wave'],refl_dict['flux'],fill_value="extrapolate")
    # mmt
    #wave_norm_mmt=4600
    #refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'juno'+'/3_Juno_27Aug09_refscl.csv', 'csv', 0, 1,
    #                      sigma_index=None, flag_index=None,
    #                      ext=1, sep=',', skiprows=0)
    #refl_mmt = normRefl(refl_mmt, wave_norm=wave_norm_mmt)
    #plt.plot(refl_mmt['wave'],refl_mmt['flux']*(f(wave_norm_mmt)/f(3800)),'k:',label='MMT',linewidth=0.5,alpha=0.5)
    plt.xlim(1680,6800)
    plt.ylim(-0.2,5)
    plt.axvline(x=3850,c='gray',lw=1,alpha=0.2)
    plt.axvline(x=4260,c='gray',lw=1,alpha=0.2)
    plt.fill_between([5000,6800], [-0.5,-0.5], [5,5], color='k', alpha=0.1, hatch='/')
    plt.ylabel('Normalized Flux')
    plt.xlabel('Wavelength')
    plt.legend()
    plt.show()

#compareAnalog()
def comparePara():
    # for 00091702002
    files = [{'name':'00091706004.pha', 'label':'Swift: RA/DEC_PNT', 'c':'r'},]
             #{'name':'00091706004_radec.pha', 'label':'Swift: by eye', 'c':'b'}]#{'name':'', 'label':,}
    for f in files:
        path = '/Users/zexixing/Research/swiftASTER/docs/'+f['name']
        spec_dict = readinSpec(path, 
                               'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                               flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
        spec_dict = normRefl(spec_dict)
        plotSepc(spec_dict['wave'], spec_dict['flux'], color=f['c'], linestyle='--',alpha=0.7, label=f['label'],
                 sigma=None, flag=None)  
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/sunSpec/sun_ref_colina96.asc.txt', 'csv', 0, 1,
                           sigma_index=None, flag_index=None,
                           ext=1, sep='    ', skiprows=0)
    spec_dict = normRefl(spec_dict)
    plotSepc(spec_dict['wave'], spec_dict['flux'], color='k', linestyle='-',alpha=1, label='colina96',
             sigma=None, flag=None)   
    plt.show()

def compareAster(aster):
    coord = '_2order_4width'
    if aster == 'juno':
        t_uvot('00091026003'+coord+'.pha',aster,'pink',linestyle='--',alpha=1,legend='00091026003')
        t_uvot('00091203002'+coord+'.pha',aster,'r',linestyle='--',alpha=1,legend='00091203002')
        t_uvot('00091204002'+coord+'.pha',aster,'y',linestyle='--',alpha=1,legend='00091204002')
        t_uvot('00091205002'+coord+'.pha',aster,'g',linestyle='--',alpha=1,legend='00091205002')
        t_uvot('00091206002'+coord+'.pha',aster,'cyan',linestyle='--',alpha=1,legend='00091206002')
        t_uvot('00091207002'+coord+'.pha',aster,'b',linestyle='--',alpha=1,legend='00091207002')
        mean = coadd_(['00091026003'+coord+'.pha','00091203002'+coord+'.pha',
                    '00091204002'+coord+'.pha','00091205002'+coord+'.pha',
                    '00091206002'+coord+'.pha'],aster)
        plotSepc(mean['wave'], np.array(mean['flux'])*1e12, color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        t_iue(['LWR01896LL.FITS','LWR05678LL.FITS','LWR05679LL.FITS','LWR05690LL.FITS','LWR06487LL.FITS'])
        plt.ylim(-0.05,1.05)
    elif aster == 'hebe':
        #['00091274001','00091510001','00091512001','00091514001','00091516001']
        t_uvot('00091274001'+coord+'.pha',aster,'pink',linestyle='--',alpha=1,legend='00091274001')
        t_uvot('00091510001'+coord+'.pha',aster,'r',linestyle='--',alpha=1,legend='00091510001')
        t_uvot('00091512001'+coord+'.pha',aster,'y',linestyle='--',alpha=1,legend='00091512001')
        t_uvot('00091514001'+coord+'.pha',aster,'g',linestyle='--',alpha=1,legend='00091514001')
        t_uvot('00091516001'+coord+'.pha',aster,'cyan',linestyle='--',alpha=1,legend='00091516001')
        mean = coadd_(['00091274001'+coord+'.pha','00091510001'+coord+'.pha',
                    '00091512001'+coord+'.pha','00091514001'+coord+'.pha',
                    '00091516001'+coord+'.pha'],aster)
        plotSepc(mean['wave'], mean['flux'], color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        plt.ylim(-0.05e-12,1.05e-12)
    elif aster == 'ceres':
        #['00091242001','00091244001','00091246001','00091248001']
        t_uvot('00091242001'+coord+'.pha',aster,'pink',linestyle='--',alpha=1,legend='00091242001')
        t_uvot('00091244001'+coord+'.pha',aster,'r',linestyle='--',alpha=1,legend='00091244001')
        t_uvot('00091246001'+coord+'.pha',aster,'y',linestyle='--',alpha=1,legend='00091246001')
        t_uvot('00091248001'+coord+'.pha',aster,'g',linestyle='--',alpha=1,legend='00091248001')
        mean = coadd_(['00091242001'+coord+'.pha','00091244001'+coord+'.pha',
                    '00091246001'+coord+'.pha','00091248001'+coord+'.pha'],aster)
        plotSepc(mean['wave'], mean['flux'], color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        plt.ylim(-0.05e-12,5e-12)
    elif aster == 'iris':
        #['00091519001','00091521001','00091523001','00091525001','00091527001']
        t_uvot('00091519001'+coord+'.pha',aster,'pink',linestyle='--',alpha=1,legend='00091519001')
        t_uvot('00091521001'+coord+'.pha',aster,'r',linestyle='--',alpha=1,legend='00091521001')
        t_uvot('00091523001'+coord+'.pha',aster,'y',linestyle='--',alpha=1,legend='00091523001')
        t_uvot('00091525001'+coord+'.pha',aster,'g',linestyle='--',alpha=1,legend='00091525001')
        t_uvot('00091527001'+coord+'.pha',aster,'cyan',linestyle='--',alpha=1,legend='00091527001')
        mean = coadd_(['00091519001'+coord+'.pha','00091521001'+coord+'.pha',
                      '00091523001'+coord+'.pha','00091525001'+coord+'.pha',
                      '00091527001'+coord+'.pha'],aster)
        plotSepc(mean['wave'], mean['flux'], color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        plt.ylim(-0.05e-12,1.05e-12)
    elif aster == 'massalia':
        #['00091541001','00091543001','00091545001','00091547001','00091549001']
        t_uvot('00091541001'+coord+'.pha',aster,'pink',linestyle='--',alpha=1,legend='00091541001')
        t_uvot('00091543001'+coord+'.pha',aster,'r',linestyle='--',alpha=1,legend='00091543001')
        t_uvot('00091545001'+coord+'.pha',aster,'y',linestyle='--',alpha=1,legend='00091545001')
        t_uvot('00091547001'+coord+'.pha',aster,'g',linestyle='--',alpha=1,legend='00091547001')
        t_uvot('00091549001'+coord+'.pha',aster,'cyan',linestyle='--',alpha=1,legend='00091549001')
        mean = coadd_(['00091541001'+coord+'.pha','00091543001'+coord+'.pha',
                      '00091545001'+coord+'.pha','00091547001'+coord+'.pha',
                      '00091549001'+coord+'.pha'],aster)
        plotSepc(mean['wave'], mean['flux'], color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        plt.ylim(-0.05e-12,1.05e-12)
    else:
        idlist = ugrismlist(aster)
        clist = np.linspace(0,1,len(idlist))
        for i in len(idlist):
            obsid = idlist[i]
            color = colorFader(clist[i])
            t_uvot(obsid+coord+'.pha',aster,color,linestyle='--',alpha=1,legend=obsid)
        coaddlist = np.array(idlist)+coord+'.pha'
        mean = coadd_(coaddlist,aster)
        plotSepc(mean['wave'], mean['flux'], color='k', label='mean',linestyle='--',linewidth=1.3,sigma=mean['sigma'])
        plt.ylim(-0.05e-12,1.05e-12)
    plt.xlim(1600,7000)
    plt.title(aster)
    plt.show()
    

def plot_refl():
    # 0009102600
    coord='_swift'
    '''
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    refl = getRefl('00091026003.pha',solar_dict)
    plt.plot(refl[0],refl[1],'r--',label='Swift: 00091026003/solar analog',linewidth=0.7)
    solar_dict = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    refl = getRefl('00091026003.pha', solar_dict)
    plt.plot(refl[0],refl[1],'r:',label='Swift: 00091026003/colina96',linewidth=0.7)
    '''
    # mean
    mean = coadd_(['00091026003'+coord+'.pha','00091203002'+coord+'.pha',
                  '00091204002'+coord+'.pha','00091205002'+coord+'.pha',
                  '00091206002'+coord+'.pha'],'juno')
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    refl = getRefl(mean, solar_dict)
    refl = redden(refl[0], refl[1], 0.5, l0=5500, l1=3000)
    plt.plot(refl[0],refl[1],'k-',label='Juno',linewidth=0.7)
    #solar_dict = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    #refl = getRefl(mean, solar_dict)
    #plt.plot(refl[0],refl[1],'k:',label='Swift: mean/colina96',linewidth=0.7)
    #solar_dict = getSolarFlux('sun_e490.csv', if_smooth=True)
    #refl = getRefl(mean, solar_dict)
    #plt.plot(refl[0],refl[1],'r--',label='Swift: mean/E490',linewidth=0.7)
    #solar_dict = getSolarFlux('sun_wehrli.csv', if_smooth=True)
    #refl = getRefl(mean, solar_dict)
    #plt.plot(refl[0],refl[1],'r-.',label='Swift: mean/Wehrli',linewidth=0.7)
    f = interpolate.interp1d(refl[0],refl[1],fill_value="extrapolate")
    #----------
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    mean = coadd_(['00091274001'+coord+'.pha','00091510001'+coord+'.pha',
                  '00091512001'+coord+'.pha','00091514001'+coord+'.pha',
                  '00091516001'+coord+'.pha'],'hebe')
    refl = getRefl(mean, solar_dict)
    plt.plot(refl[0],refl[1],'r-',label='Hebe',linewidth=0.7)
    mean = coadd_(['00091242001'+coord+'.pha','00091244001'+coord+'.pha',
                  '00091246001'+coord+'.pha','00091248001'+coord+'.pha'],'ceres')
    refl = getRefl(mean, solar_dict)
    plt.plot(refl[0],refl[1],'b-',label='Ceres',linewidth=0.7)
    mean = coadd_(['00091519001'+coord+'.pha','00091521001'+coord+'.pha',
                  '00091523001'+coord+'.pha','00091525001'+coord+'.pha',
                  '00091527001'+coord+'.pha'],'iris')
    refl = getRefl(mean, solar_dict)
    plt.plot(refl[0],refl[1],'g-',label='Iris',linewidth=0.7)
    mean = coadd_(['00091541001'+coord+'.pha','00091543001'+coord+'.pha',
                  '00091545001'+coord+'.pha','00091547001'+coord+'.pha',
                  '00091549001'+coord+'.pha'],'massalia')
    refl = getRefl(mean, solar_dict)
    plt.plot(refl[0],refl[1],'y-',label='Massalia',linewidth=0.7)
    #----------
    '''
    # iue
    refl_iue = readinSpec('/Users/zexixing/Dropbox/PDART Small Bodies UV Archive/IUE/data/reflectance/3 Juno/juno_avg_refl.tab', 'txt', 0, 1,
                      sigma_index=None, flag_index=None,
                      ext=1, sep=' ', skiprows=0)
    refl_iue = normRefl({'wave':np.array(refl_iue['wave'])*10,'flux':np.array(refl_iue['flux'])}, wave_norm=3000)
    plt.plot(refl_iue['wave'],refl_iue['flux']*(f(3000)/f(5500)),'r-',label='IUE',linewidth=0.5,alpha=0.5)
    
    # mmt
    refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/3_Juno_27Aug09_refscl.csv', 'csv', 0, 1,
                          sigma_index=None, flag_index=None,
                          ext=1, sep=',', skiprows=0)
    refl_mmt = normRefl(refl_mmt)
    plt.plot(refl_mmt['wave'],refl_mmt['flux'],'b-',label='MMT',linewidth=0.5,alpha=0.5)
    
    # smass
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass2.txt', 
                            'txt', 0, 1,  sigma_index=2,ext=1, sep='\t', skiprows=0)
    spec_dict = binSpec(np.array(spec_dict['wave'])*1e4, np.array(spec_dict['flux']), binGenerate(np.array(spec_dict['wave'])*1e4, 10))
    refl = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=5500)
    plotSepc(spec_dict['wave'], np.array(spec_dict['flux']), color='g', label='smass2',linestyle='--')
    '''
    '''
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass1.txt', 
                            'txt', 0, 1,  sigma_index=None,ext=1, sep=' ', skiprows=0)
    spec_dict = binSpec(spec_dict['wave'], np.array(spec_dict['flux']), binGenerate(spec_dict['wave'], 10), sigma=None)#np.array(spec_dict['sigma'])*norm)
    refl = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=5500)
    plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(5500)/f(3000)), color='g', label='smass1',linestyle=':',sigma=None)#spec_dict['sigma'])

    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass3.txt', 
                            'txt', 0, 1,  sigma_index=2,ext=1, sep='     ', skiprows=0)
    spec_dict = binSpec(np.array(spec_dict['wave'])*1e4, np.array(spec_dict['flux']), binGenerate(np.array(spec_dict['wave'])*1e4, 10))#, sigma=np.array(spec_dict['sigma'])*norm)
    refl = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=9000)
    plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(9000)/f(3000)), color='k', label='smass3',linestyle='-.')#,sigma=spec_dict['sigma'])
    '''
    plt.xlim(1600,7000)
    plt.ylim(-0.5,1.5)
    plt.ylabel('Normalized reflectance')
    plt.xlabel('wavelength (Angstrom)')
    plt.legend()
    plt.show()

def plot_refl_mean(swift=True, mmt=False, iue=False, 
                   juno=False, hebe=False, ceres=False, iris=False, massalia=False, allaster=True):
    coord='_smeargauss'
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    #solar_dict = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    if swift == True:
        if juno == True or allaster == True:
            refl = coadd_([getRefl('00091026003'+coord+'.pha', solar_dict, 'juno'),
                        getRefl('00091203002'+coord+'.pha', solar_dict, 'juno'),
                        getRefl('00091204002'+coord+'.pha', solar_dict, 'juno'),
                        getRefl('00091205002'+coord+'.pha', solar_dict, 'juno'),
                        getRefl('00091206002'+coord+'.pha', solar_dict, 'juno')])
            #refl = redden(refl['wave'], refl['flux'], 0.5, l0=5500, l1=3000)
            #wave, flux = refl
            refl = normRefl(refl)
            wave, flux = refl['wave'], refl['flux']
            plt.plot(wave,flux,'k-',label='Juno',linewidth=0.7)
            f = interpolate.interp1d(wave,flux,fill_value="extrapolate")
        if hebe == True or allaster == True: 
            refl = coadd_([getRefl('00091274001'+coord+'.pha', solar_dict, 'hebe'),
                        getRefl('00091510001'+coord+'.pha', solar_dict, 'hebe'),
                        getRefl('00091512001'+coord+'.pha', solar_dict, 'hebe'),
                        getRefl('00091514001'+coord+'.pha', solar_dict, 'hebe'),
                        getRefl('00091516001'+coord+'.pha', solar_dict, 'hebe')])
            plt.plot(refl['wave'],refl['flux'],'r-',label='Hebe',linewidth=0.7)
            f_hebe = interpolate.interp1d(refl['wave'],refl['flux'],fill_value="extrapolate")
        if ceres == True or allaster == True: 
            refl = coadd_([getRefl('00091242001'+coord+'.pha', solar_dict, 'ceres'),
                        getRefl('00091244001'+coord+'.pha', solar_dict, 'ceres'),
                        getRefl('00091246001'+coord+'.pha', solar_dict, 'ceres'),
                        getRefl('00091248001'+coord+'.pha', solar_dict, 'ceres')])
            plt.plot(refl['wave'],refl['flux'],'b-',label='Ceres',linewidth=0.7)
            f_ceres = interpolate.interp1d(refl['wave'],refl['flux'],fill_value="extrapolate")
        if iris == True or allaster == True: 
            refl = coadd_([getRefl('00091519001'+coord+'.pha', solar_dict, 'iris'),
                        getRefl('00091521001'+coord+'.pha', solar_dict, 'iris'),
                        getRefl('00091523001'+coord+'.pha', solar_dict, 'iris'),
                        getRefl('00091525001'+coord+'.pha', solar_dict, 'iris'),
                        getRefl('00091527001'+coord+'.pha', solar_dict, 'iris')])
            plt.plot(refl['wave'],refl['flux'],'g-',label='Iris',linewidth=0.7)
            f_iris = interpolate.interp1d(refl['wave'],refl['flux'],fill_value="extrapolate")
        if massalia == True or allaster == True:
            refl = coadd_([getRefl('00091541001'+coord+'.pha', solar_dict, 'massalia'),
                        getRefl('00091543001'+coord+'.pha', solar_dict, 'massalia'),
                        getRefl('00091545001'+coord+'.pha', solar_dict, 'massalia'),
                        getRefl('00091547001'+coord+'.pha', solar_dict, 'massalia'),
                        getRefl('00091549001'+coord+'.pha', solar_dict, 'massalia')])
            plt.plot(refl['wave'],refl['flux'],'y-',label='Massalia',linewidth=0.7)
            f_massalia = interpolate.interp1d(refl['wave'],refl['flux'],fill_value="extrapolate")
    
    # mmt
    if mmt == True:
        if juno == True or allaster == True:
            refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'juno'+'/3_Juno_27Aug09_refscl.csv', 'csv', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=',', skiprows=0)
            refl_mmt = normRefl(refl_mmt, wave_norm=3800)
            plt.plot(refl_mmt['wave'],refl_mmt['flux']*(f(3800)/f(3800)),'r--',label='MMT',linewidth=0.5,alpha=0.5)
        if hebe == True or allaster == True:
            refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'hebe'+'/6_Hebe_18Jun13_refscl.csv', 'csv', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=',', skiprows=0)
            refl_mmt = normRefl(refl_mmt)
            plt.plot(refl_mmt['wave'],refl_mmt['flux'],'r--',label='MMT Hebe',linewidth=0.5,alpha=0.5)
        if iris == True or allaster == True:
            refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'iris'+'/7_Iris_11Dec06_refscl.csv', 'csv', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=',', skiprows=0)
            refl_mmt = normRefl(refl_mmt)
            plt.plot(refl_mmt['wave'],refl_mmt['flux'],'g--',label='MMT Iris',linewidth=0.5,alpha=0.5)    
        if massalia == True or allaster == True:
            refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'massalia'+'/20_Massalia_UT11Dec06_refscl.csv', 'csv', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=',', skiprows=0)
            refl_mmt = normRefl(refl_mmt)
            plt.plot(refl_mmt['wave'],refl_mmt['flux'],'y--',label='MMT Massalia',linewidth=0.5,alpha=0.5)   

    # iue
    if iue == True:
        if juno or allaster == True:
            refl_iue = readinSpec('/Users/zexixing/Dropbox/PDART Small Bodies UV Archive/IUE/data/reflectance/3 Juno/juno_avg_refl.tab', 'txt', 0, 1,
                            sigma_index=None, flag_index=None,
                            ext=1, sep=' ', skiprows=0)
            refl_iue = normRefl({'wave':np.array(refl_iue['wave'])*10,'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f(2600)/f(3800)),'b:',label='IUE Juno',linewidth=0.5,alpha=0.5)
        if hebe == True or allaster == True:
            refl_iue = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'hebe'+'/hebe_iue.txt', 'txt', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep='', skiprows=0)
            refl_iue = normRefl({'wave':np.array(refl_iue['wave']),'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f_hebe(2600)/f_hebe(5500)),'r:',label='IUE Hebe',linewidth=0.5,alpha=0.5)
        if iris == True or allaster == True:
            refl_iue = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'iris'+'/iris_iue.txt', 'txt', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep='', skiprows=0)
            refl_iue = normRefl({'wave':np.array(refl_iue['wave']),'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f_iris(2600)/f_iris(5500)),'g:',label='IUE Iris',linewidth=0.5,alpha=0.5)
        if massalia == True or allaster == True:
            refl_iue = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'massalia'+'/massalia_iue.txt', 'txt', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep='', skiprows=0)
            refl_iue = normRefl({'wave':np.array(refl_iue['wave']),'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f_massalia(2600)/f_massalia(5500)),'y:',label='IUE Massalia',linewidth=0.5,alpha=0.5)
        if ceres == True or allaster == True:
            refl_iue = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'ceres'+'/ceres_iue.txt', 'txt', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep='', skiprows=0)
            refl_iue = normRefl({'wave':np.array(refl_iue['wave']),'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f_ceres(2600)/f_ceres(5500)),'b:',label='IUE Ceres',linewidth=0.5,alpha=0.5)
    # smass - temp
    if juno==True:
        f_mmt = interpolate.interp1d(refl_mmt['wave'],refl_mmt['flux'],fill_value="extrapolate")
        aster = 'juno'
        spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass2.txt', 
                                'txt', 0, 1,  sigma_index=2,ext=1, sep='\t', skiprows=0)
        spec_dict = binSpec(np.array(spec_dict['wave'])*1e4, np.array(spec_dict['flux']), binGenerate(np.array(spec_dict['wave'])*1e4, 10))
        spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=5500)
        plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(5500)/f(3800)), color='g', label='smass2',linestyle='--')

        spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass1.txt', 
                                'txt', 0, 1,  sigma_index=None,ext=1, sep=' ', skiprows=0)
        spec_dict = binSpec(spec_dict['wave'], np.array(spec_dict['flux']), binGenerate(spec_dict['wave'], 10), sigma=None)#np.array(spec_dict['sigma'])*norm)
        spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=5500)
        plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(5500)/f(3800)), color='g', label='smass1',linestyle=':',sigma=None)#spec_dict['sigma'])
    plt.xlim(1680,6800)
    plt.ylim(-0.5,2)
    plt.ylabel('Normalized reflectance')
    plt.xlabel('wavelength (Angstrom)')
    #plt.axvline(x=6000,c='gray',lw=155,alpha=0.2)
    plt.fill_between([5000,6800], [-0.5,-0.5], [2,2], color='k', alpha=0.3, hatch='/')
    plt.legend()
    plt.show()


def plot_refl_sun():
    coord='_swift'
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=False)
    refl = getRefl('00091026003'+coord+'.pha', solar_dict, 'juno')
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'r--',label='Juno 00091026003',linewidth=0.5,alpha=0.5)
    refl = getRefl('00091203002'+coord+'.pha', solar_dict, 'juno')
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'b--',label='Juno 00091203002',linewidth=0.5,alpha=0.5)
    refl = getRefl('00091204002'+coord+'.pha', solar_dict, 'juno')
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'y--',label='Juno 00091204002',linewidth=0.5,alpha=0.5)
    refl = getRefl('00091205002'+coord+'.pha', solar_dict, 'juno')
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'g--',label='Juno 00091205002',linewidth=0.5,alpha=0.5)
    refl = getRefl('00091206002'+coord+'.pha', solar_dict, 'juno')
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'c--',label='Juno 00091206002',linewidth=0.5,alpha=0.5) 
    solar_path = '/Users/zexixing/Research/swiftASTER/docs/'+'00091706004.pha'
    solar_spec_dict = readinSpec(solar_path, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    solar_spec_dict=coadd_([solar_spec_dict])
    solar_spec_dict = normRefl(solar_spec_dict)
    plt.plot(solar_spec_dict['wave'],solar_spec_dict['flux'],'k:',label='solar analog (smooth)',linewidth=0.8)
    solar_spec_dict = readinSpec(solar_path, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    solar_spec_dict = normRefl(solar_spec_dict)
    plt.plot(solar_spec_dict['wave'],solar_spec_dict['flux'],'k-',label='solar analog',linewidth=0.8)
    plt.xlim(2150,6800)
    plt.ylim(-0.5,1.5)
    plt.ylabel('Normalized reflectance')
    plt.xlabel('wavelength (Angstrom)')
    plt.axvline(x=6000,c='gray',lw=155,alpha=0.2)
    plt.legend()
    plt.show()

def refl_smooth():
    coord='_swift'
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    # 132
    refl = coadd_([getRefl('00091274001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091510001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091512001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091514001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091516001'+coord+'.pha', solar_dict, 'hebe')]) 
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'r-',label='smooth, divide, average',linewidth=0.7)
    # 123
    mean = coadd_(['00091274001'+coord+'.pha','00091510001'+coord+'.pha',
                  '00091512001'+coord+'.pha','00091514001'+coord+'.pha',
                  '00091516001'+coord+'.pha'],'hebe')
    refl = getRefl(mean, solar_dict)
    plt.plot(refl['wave'],refl['flux'],'g-',label='smooth, average, divide',linewidth=0.7)
    # 231
    solar_dict = getSolarFlux('00091706004.pha', if_smooth=False)
    mean = coadd_(['00091274001'+coord+'.pha','00091510001'+coord+'.pha',
                  '00091512001'+coord+'.pha','00091514001'+coord+'.pha',
                  '00091516001'+coord+'.pha'],'hebe')
    refl = getRefl(mean, solar_dict)
    refl = coadd_([refl])
    plt.plot(refl['wave'],refl['flux'],'b-',label='average, divide, smooth',linewidth=0.7)
    # 3(21)
    refl = coadd_([getRefl('00091274001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091510001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091512001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091514001'+coord+'.pha', solar_dict, 'hebe'),
                getRefl('00091516001'+coord+'.pha', solar_dict, 'hebe')])
    plt.plot(refl['wave'],refl['flux'],'y-',label='divide, smooth+average',linewidth=0.7)
    solar_path = '/Users/zexixing/Research/swiftASTER/docs/'+'00091706004.pha'
    solar_spec_dict = readinSpec(solar_path, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    solar_spec_dict=coadd_([solar_spec_dict])
    solar_spec_dict = normRefl(solar_spec_dict)
    plt.plot(solar_spec_dict['wave'],0.7*solar_spec_dict['flux'],'k:',label='solar analog (smooth)',linewidth=0.8)
    solar_spec_dict = readinSpec(solar_path, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    solar_spec_dict = normRefl(solar_spec_dict)
    plt.plot(solar_spec_dict['wave'],0.7*solar_spec_dict['flux'],'k-',label='solar analog',linewidth=0.8)
    plt.xlim(2150,6800)
    plt.ylim(-0.5,1.5)
    plt.ylabel('Normalized reflectance')
    plt.xlabel('wavelength (Angstrom)')
    plt.axvline(x=6000,c='gray',lw=155,alpha=0.2)
    plt.legend()
    plt.show()

def asterid(aster):
    aster_dict = {'ceres':1, 'juno':3, 'hebe':6, 'dembowska':349,
                  'eros':433, 'eunomia':15, 'flora':8, 'hygiea':10,
                  'iris':7, 'lutetia':21, 'massalia':20, 'nysa':44,
                  'pallas':2, 'themis':24, 'toutatis':4179, 'vesta':4,
                  'scheila':596, 'yu':55}
    return aster_dict[aster]

def swift16flux(asterlist,remove2nd=False,iue=True,part=False):
    c_dict = {'ceres':'firebrick', 'pallas':'red', 'juno':'darksalmon', 'vesta':'gold', 'hebe':'olivedrab', 'iris':'mediumspringgreen', 'flora':'lightseagreen', 
              'hygiea':'deepskyblue', 'eunomia':'royalblue', 'massalia':'navy', 'lutetia':'blue', 'themis':'mediumpurple', 'nysa':'darkorchid',
              'dembowska':'plum', 'eros':'mediumvioletred', 'scheila':'palevioletred', 'toutatis':'sandybrown', 
             }
    c_pool = ['firebrick','navy','blue','mediumspringgreen','darkorchid','gold','lightseagreen','red','olivedrab','sandybrown',
              'mediumpurple','deepskyblue','royalblue','mediumvioletred','palevioletred','plum','darksalmon']
    import plotly as py
    pympl = py.offline.plot_mpl
    plt.figure(figsize=(6,4))
    coord = '_radec'
    clist = np.linspace(0,1,len(asterlist))
    remove_list = ['00091503001', #flora -- star
                   '00091223001','00091225001','00091229001', #pallas -- star(compact)
                   '00091198002', #'00091022002','00091197002', vesta -- star
                   '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                   '00091538001','00091540001',#nysa -- star(compact)
                   '00091591001', #themis -- wrong position
                   '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                   ]
    for i in range(0,len(asterlist)):
        aster = asterlist[i]
        color = c_dict[aster] #colorFader(mix=clist[i])
        idlist = ugrismlist(aster)
        for obsid_rm in remove_list:
            try:    idlist.remove(obsid_rm)
            except:    pass
        if len(asterlist)==1:
            label='mean'
            plt.title(str(int(asterid(aster)))+' '+aster)
            for j in range(0,len(idlist)):
                obsid = idlist[j]
                color = c_pool[j]
                clist = np.linspace(0,1,len(idlist))
                color = colorFader(mix=clist[j],c1='green',c2='blue')
                #color = 'lightgray'
                if part == False:
                    t_uvot(obsid+coord+'.pha',aster,color,linestyle='-',alpha=0.5,remove2nd=remove2nd)#,legend=obsid)
                else:
                    if obsid in part:
                        t_uvot(obsid+coord+'.pha',aster,color,linestyle='-',alpha=0.7,legend=obsid,remove2nd=remove2nd)
            color = 'k'
        else:
            label=aster
        coaddlist = [obsid+coord+'.pha' for obsid in idlist]
        mean = coadd(coaddlist,aster,remove2nd=remove2nd)
        plotSepc(mean['wave'], np.array(mean['flux'])*1e12, color=color, linestyle='-', linewidth=1, label=label)
        rela_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
        if iue:
            if (aster!='themis') and (aster!='scheila') and (aster!='toutatis'):
                if aster == 'juno':
                    iue = t_iue(['LWR01896LL.FITS','LWR05678LL.FITS','LWR05679LL.FITS','LWR05690LL.FITS','LWR06487LL.FITS'],plot=False)
                else:
                    data_path_iue = rela_path + aster+'_raw_for_zexi.txt'
                    iue = readinSpec(data_path_iue, 'txt', 0, 1, sep='')
                wave_iue = np.array(iue['wave'])
                flux_iue = np.array(iue['flux'])
                if len(asterlist)==1: label = 'iue'
                else: label = None
                plt.plot(wave_iue[wave_iue<3200],flux_iue[wave_iue<3200]*1e12,
                        color=color, linestyle='--', linewidth=1, label=label)
    plt.ylim(-0.05,1.05) #--TODO:
    plt.xlim(3200,5250)
    plt.ylim(0.1,1.2)
    #plt.yscale('log')
    #plt.xscale('log')
    fig = plt.gcf()
    plt.ylabel('flux (10^-12 erg s-1 cm-2 A-1)')
    plt.xlabel('wavelength (A)')
    #plot_url = pympl(fig,filename=r'output_plot.html',show_link=False,resize=True)
    '''
    # sun
    solar_path = '/Users/zexixing/Research/swiftASTER/docs/solar_analog/'+'00091706004_2order.pha'
    spec_dict = readinSpec(solar_path, 
                            'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                            flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    f = interpolate.interp1d(mean['wave'], np.array(mean['flux'])*1e12)
    spec_dict = normRefl(spec_dict)
    plotSepc(spec_dict['wave'], spec_dict['flux']*f(3900), color='k', linestyle=':',alpha=1, label='solar analog',#analog_list[i][2:13],
                sigma=None, flag=None)
    '''
    plt.legend()
    plt.show()

def rot(image, xy, angle):
    from scipy.ndimage import rotate
    im_rot = rotate(image,angle,reshape=False,order=1) 
    org_center = (np.array(image.shape[:2][::-1])-1)/2.
    rot_center = (np.array(im_rot.shape[:2][::-1])-1)/2.
    org = xy-org_center
    a = np.deg2rad(angle)
    new = np.array([org[0]*np.cos(a) + org[1]*np.sin(a),
            -org[0]*np.sin(a) + org[1]*np.cos(a) ])
    return im_rot, new+rot_center

def asterMotion(aster_dict,output_name,coord):
    # for all asteroids
    aster = {"w": "ASTEROID"}
    obs = {"w": "OBS_ID"}
    start_t = {"r": "DATE-OBS"}
    end_t = {"r": "DATE-END"}
    mid_t = {"w": "MIDTIME"}
    exp = {"w": "EXPOSURE",
           "r": "EXPOSURE"}
    angle = {"w": "ANGLE",
             "r": "ANGLE"}
    m_para = {"w": "MOTION_P"}
    m_vert = {"w": "MOTION_V"}
    m_all = {"w": "MOTION_A"}
    track = {"w": "TRACKWIDTH",
             "r1": "WIDTH_UP",
             "r2": "WIDTH_DW"}
    sigma = {"w": "SIGMA"}
    para_dict = {'aster':aster,
                 'obs':obs,
                 'start_t':start_t,
                 'end_t':end_t,
                 'mid_t':mid_t,
                 'exp':exp,
                 'angle':angle,
                 'm_all':m_all,
                 'm_para':m_para,
                 'm_vert':m_vert,
                 'track':track,
                 'sigma':sigma}
    output_path = get_path('../docs/'+output_name)
    f = open(output_path, 'a')
    #f.write(  aster["w"] + ' '
    #        + obs["w"] + ' '
    #        + mid_t["w"] + ' '
    #        + exp["w"] + ' '
    #        + angle["w"] + ' '
    #        + m_para["w"] + ' '
    #        + m_vert["w"] + ' '
    #        + sigma["w"] + ' '
    #        + track["w"] + '\n')
    for aster_name in aster_dict.keys():
        asterMotion_single(aster_name,aster_dict[aster_name],para_dict,coord,f)
    f.close()

def asterMotion_single(aster_name,aster_id,para_dict,coord,f):
    from astropy.time import Time
    from astropy.wcs import WCS
    aster = para_dict['aster']
    obs = para_dict['obs']
    start_t = para_dict['start_t']
    end_t = para_dict['end_t']
    mid_t = para_dict['mid_t']
    exp = para_dict['exp']
    angle = para_dict['angle']
    m_all = para_dict['m_all']
    m_para = para_dict['m_para']
    m_vert = para_dict['m_vert']
    track = para_dict['track']
    sigma = para_dict['sigma']
    aster["v"] = aster_name
    # obtaion obsid list
    obsid_list = ugrismlist(aster_name,remove=['00013947002', '00013947003'],filt='uvgrism')
    print(obsid_list)
    # loop for every obs
    for obsid in obsid_list:
        #if obsid=='00091026003':
        #    pass
        #else:
        #    break
        print(obsid)
        obs["v"] = obsid
        obs_name = obsid+coord+'.pha'
        obs_path = get_path('../docs/'+aster_name+'/'+obs_name)
        obs_hdul = fits.open(obs_path)
        angle["v"] = obs_hdul[0].header[angle["r"]]
        exp["v"] = obs_hdul[1].header[exp["r"]]
        start_t["v"] = Time(obs_hdul[1].header[start_t["r"]])
        end_t["v"] = Time(obs_hdul[1].header[end_t["r"]])
        dt = end_t["v"] - start_t["v"]
        mid_t["v"] = start_t["v"] + 1/2*dt
        track["v"] = np.mean(obs_hdul[2].data['WIDTH_UP']-obs_hdul[2].data['WIDTH_DW'])
        sigma["v"] = track["v"]/8
        obs_name = 'sw'+obsid+'ugu_dt.img'
        obs_path = get_path('../data/'+aster_name+'/'+obsid+'/uvot/image/'+obs_name)
        obs_hdul = fits.open(obs_path)
        obs_img = obs_hdul[1].data
        w = WCS(obs_hdul[1].header,key='S',relax=True)
        obj_start = Horizons(id=aster_id,
                             location='@swift', 
                             epochs=start_t["v"].jd)
        eph_start = obj_start.ephemerides()[0]
        ra_start = eph_start['RA']
        dec_start = eph_start['DEC']
        px_start,py_start = w.wcs_world2pix(ra_start,dec_start,0)
        obj_end = Horizons(id=aster_id,
                           location='@swift', 
                           epochs=end_t["v"].jd)
        eph_end = obj_end.ephemerides()[0]
        ra_end = eph_end['RA']
        dec_end = eph_end['DEC']
        px_end,py_end = w.wcs_world2pix(ra_end,dec_end,0)
        # dist
        dx = px_end - px_start
        dy = py_end - py_start
        if angle["v"]<40:
            theta = (180-angle["v"])*np.pi/180
        else:
            theta = (angle["v"])*np.pi/180
        c1 = dx/np.cos(theta)
        c2 = dy - dx*np.tan(theta)
        m_all["v"] = np.sqrt((px_start-px_end)**2+(py_start-py_end)**2)
        m_para["v"] = c1 + c2*np.sin(theta)
        m_vert["v"] = c2*np.cos(theta)
        # another method
        #---new,new_start = rot(obs_img, (px_start-1,py_start-1), 180-angle["v"])
        #---new,new_end = rot(obs_img, (px_end-1,py_end-1), 180-angle["v"])
        #---new_para = new_end[0]-new_start[0]
        #---new_vert = new_end[1]-new_start[1]
        # write
        f.write(  aster["v"] + ' '
                + obs["v"] + ' '
                + str(mid_t["v"]) + ' '
                + str(exp["v"]) + ' '
                + str(angle["v"]) + ' '
                + str(m_para["v"]) + ' '
                + str(m_vert["v"]) + ' '
                + '%.4f'%sigma["v"] + ' '
                + '%.4f'%track["v"]+ '\n')
        
def motion_plot(plotly):
    import plotly as py
    pympl = py.offline.plot_mpl
    fig = plt.figure(figsize=(6,5))
    file_name = 'motion_log.txt'
    file_path = get_path('../docs/'+file_name)
    #f = open(file_path)
    #for line in f.readlines()[1:]:
    data = pd.read_csv(file_path, sep=' ', header=0)
    data_group = data.groupby('ASTEROID')
    c_dict = {'ceres':'firebrick', 'pallas':'red', 'juno':'darksalmon', 'vesta':'gold', 'hebe':'olivedrab', 'iris':'mediumspringgreen', 'flora':'lightseagreen', 
              'hygiea':'deepskyblue', 'eunomia':'royalblue', 'massalia':'navy', 'lutetia':'blue', 'themis':'mediumpurple', 'nysa':'darkorchid',
              'dembowska':'plum', 'eros':'mediumvioletred', 'scheila':'palevioletred', 'toutatis':'sandybrown', 
             }
    #for name, group in data_group:
    #    motion_perp = np.array(group['MOTION_V'])
    #    motion_para = np.array(group['MOTION_P'])
    #    plt.plot(motion_para, motion_perp, 'o', color=c_dict[name], label = name, MarkerSize=3)
    for i in range(0,len(data)):
        motion_para = data['MOTION_P'][i]
        motion_perp = data['MOTION_V'][i]
        aster = data['ASTEROID'][i]
        obsid = data['OBS_ID'][i]
        if plotly == True: label = aster+':'+str(int(obsid))[2:]
        else: label = aster
        plt.plot(motion_para, motion_perp, 'o', color=c_dict[aster], label = label, MarkerSize=3)
    sigma = 3.5907
    plt.hlines(sigma*8, -30, 180, ls='--')
    plt.hlines(-sigma*8, -30, 180, ls='--')
    plt.fill_between([-30,180],[sigma,sigma],[-sigma,-sigma],color='gray',alpha=0.2,hatch='/')
    plt.xlim(-30,180)
    plt.ylim(-160,55)
    plt.xlabel('parallal')
    plt.ylabel('perpendicular')
    if plotly == True: 
        plot_url = pympl(fig,filename=r'motion.html',show_link=False,resize=True)
    else: 
        #plt.legend()
        plt.show()
    
#motion_plot(plotly=True)

def average_refl(aster):
    solar_spec_dict = getSolarFlux('00091706004.pha', if_smooth=True)
    #solar_spec_dict['wave'] = solar_spec_dict['wave']-15
    solar_spec_dict = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    solar_spec_dict = coadd_([solar_spec_dict],binbox=20)
    idlist = ugrismlist(aster)
    remove_list = ['00091503001', #flora -- star
                   '00091223001','00091225001','00091229001', #pallas -- star(compact)
                   '00091198002', #'00091022002','00091197002', vesta -- star
                   '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                   '00091540001','00091538001', #nysa -- star(compact)
                   '00091591001', #themis -- wrong position
                   '00091207002', #juno -- star in slit
                   ]
    for obsid_rm in remove_list:
        try:    idlist.remove(obsid_rm)
        except:    pass
    # loop
    adjname = '_smeargauss'
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
    flux_dict = coadd_(flux_list,aster)
    #flux_dict['wave']=flux_dict['wave']-15
    refl_dict = getRefl(flux_dict, solar_spec_dict, aster=aster)
    refl = normRefl(refl_dict)
    wave, flux = refl['wave'], refl['flux']
    plt.plot(wave,flux,'k-',label='Swift (colina96)',linewidth=0.5)
    f = interpolate.interp1d(wave,flux,fill_value="extrapolate")
    # mmt
    wave_norm_mmt=4600
    refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+'juno'+'/3_Juno_27Aug09_refscl.csv', 'csv', 0, 1,
                          sigma_index=None, flag_index=None,
                          ext=1, sep=',', skiprows=0)
    refl_mmt = normRefl(refl_mmt, wave_norm=wave_norm_mmt)
    plt.plot(refl_mmt['wave'],refl_mmt['flux']*(f(wave_norm_mmt)/f(3800)),'r:',label='MMT',linewidth=0.5,alpha=0.5)
    #iue
    refl_iue = readinSpec('/Users/zexixing/Dropbox/PDART Small Bodies UV Archive/IUE/data/reflectance/3 Juno/juno_avg_refl.tab', 'txt', 0, 1,
                          sigma_index=None, flag_index=None,
                          ext=1, sep=' ', skiprows=0)
    refl_iue = normRefl({'wave':np.array(refl_iue['wave'])*10,'flux':np.array(refl_iue['flux'])}, wave_norm=2600)
    plt.plot(refl_iue['wave'],refl_iue['flux']*(f(2600)/f(3800)),'b:',label='IUE',linewidth=0.5,alpha=0.5)
    #smass
    '''
    wave_norm_smass=4600
    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass2.txt', 
                            'txt', 0, 1,  sigma_index=2,ext=1, sep='\t', skiprows=0)
    spec_dict = binSpec(np.array(spec_dict['wave'])*1e4, np.array(spec_dict['flux']), binGenerate(np.array(spec_dict['wave'])*1e4, 10))
    spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=wave_norm_smass)
    plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(wave_norm_smass)/f(3800)), color='g', label='smass2',linestyle='--')

    spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass1.txt', 
                            'txt', 0, 1,  sigma_index=None,ext=1, sep=' ', skiprows=0)
    spec_dict = binSpec(spec_dict['wave'], np.array(spec_dict['flux']), binGenerate(spec_dict['wave'], 10), sigma=None)#np.array(spec_dict['sigma'])*norm)
    spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=wave_norm_smass)
    plotSepc(spec_dict['wave'], np.array(spec_dict['flux'])*(f(wave_norm_smass)/f(3800)), color='g', label='smass1',linestyle=':',sigma=None)#spec_dict['sigma'])
    '''
    plt.xlim(1680,6800)
    plt.ylim(0,2)
    plt.axvline(x=3850,c='gray',lw=1,alpha=0.2)
    plt.axvline(x=4260,c='gray',lw=1,alpha=0.2)
    plt.ylabel('Normalized reflectance')
    plt.xlabel('wavelength (Angstrom)')
    #plt.axvline(x=6000,c='gray',lw=155,alpha=0.2)
    plt.fill_between([5000,6800], [-0.5,-0.5], [2,2], color='k', alpha=0.1, hatch='/')
    plt.legend() 
    plt.title(aster)
    plt.show()

def standard_star(star, obsid, obs_name=None, mod_name=None):
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
    # swift
    sw_name = obsid+'_default.pha'
    sw_path = docsdir+sw_name
    sw = fits.open(sw_path)
    sw_data = sw[2].data
    sw_time = (sw[2].header['TSTART']+sw[2].header['TSTOP'])/2
    fscale = (sw_time-126230400)/(12.6*365.26*86400)
    flux_sw = sw_data['FLUX']
    wave_sw = sw_data['LAMBDA']
    plt.plot(wave_sw, flux_sw, 'k-',label='Swift')
    plt.plot(wave_sw, flux_sw*fscale, 'k--', alpha=0.3, label='Swift (no fscale)')
    # others
    if obs_name != None:
        obs_path = docsdir+obs_name
        obs = fits.open(obs_path)
        obs_data = obs[1].data
        flux_obs = obs_data['FLUX']
        wave_obs = obs_data['WAVELENGTH']
        plt.plot(wave_obs, flux_obs, 'b-',label='other')
    # model
    if mod_name != None:
        mod_path = docsdir+mod_name
        mod = fits.open(mod_path)
        mod_data = mod[1].data
        flux_mod = mod_data['FLUX']
        wave_mod = mod_data['WAVELENGTH']
        plt.plot(wave_mod, flux_mod, 'r-',label='model')
    plt.xlim(800,7000)
    plt.ylim(0,20.3e-12)
    plt.title(star+' '+obsid)
    plt.legend()
    plt.show()

mod_dict = {'wd1057':['wd1057_719_stisnic_008.fits','wd1057_719_mod_006.fits'],
            'p177':['p177d_stisnic_008.fits','p177d_001.fits'],#'p177d_mod_003.fits'],
            'bd25':['bd_25d4655_002.fits'],
            'agk81':['agk_81d266_005.fits','agk_81d266_stisnic_007.fits'],
            'gd153':['gd153_fos_003.fits','gd153_mod_011.fits','gd153_stiswfcnic_002.fits'],
            'p41':['p041c_stisnic_008.fits','p041c_001.fits'],#'kurucz_kp00_7250.fits'],#'p041c_mod_003.fits'],
            'wd0320':['wd0320_539_mod_001.fits','wd0320_539_stis_005.fits'],
            'wr1':[],
            'wr4':[],
            'wr52':['iue'],
            'wr86':[],
            'wr121':[],
            'wd1657':['wd1657_343_mod_007.fits','wd1657_343_stiswfcnic_003.fits'],
           }

def read_mod(star,mod_dict):
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
    mod_list = mod_dict[star]
    if mod_list == None:
        return None
    elif mod_list[0] == 'iue':
        import glob
        mod_list = glob.glob(docsdir+'iue/'+'*.fits')
        flux_list = []
        wave_list = []
        sigma_list = []
        all_list = []
        for mod_path in mod_list:
            mod_data = fits.open(mod_path)[1].data
            flux_list.append(mod_data['FLUX'][0])
            wave_list.append(mod_data['WAVE'][0])
            sigma_list.append(mod_data['SIGMA'][0])
        spec_dict = mergeSpec(wave_list, flux_list, sigma_list=sigma_list)
        spec_dict = binSpec(spec_dict['wave'], spec_dict['flux'], binGenerate(spec_dict['wave'], 5), sigma=None, flag=None)
        mod_dict = {'wave':spec_dict['wave'],'flux':spec_dict['flux'],'type':'iue'}
        all_list.append(mod_dict)
    else:
        all_list = []
        for mod_name in mod_list:
            mod_path = docsdir+mod_name
            mod = fits.open(mod_path)
            mod_data = mod[1].data
            if 'kurucz' in mod_name:
                flux_mod = mod_data['g40']
            else:
                flux_mod = mod_data['FLUX']
            wave_mod = mod_data['WAVELENGTH']
            if 'mod' in mod_name:
                type_mod = 'model'
            elif 'stis' in mod_name:
                type_mod = 'stis'
            elif 'kurucz' in mod_name:
                type_mod = 'kurucz'
            else:
                type_mod = 'other'
            mod_dict = {'wave':wave_mod,'flux':flux_mod,'type':type_mod}
            all_list.append(mod_dict)
    return all_list
            
def star_compare(star,mod_dict):
    # mod spec
    c_dict = {'stis':'b','model':'r','other':'y','iue':'b','kurucz':'r'}
    mod_list = read_mod(star,mod_dict)
    # star spec
    import glob
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
    clock_path_list = glob.glob(docsdir+'000*_default.pha')
    normi_path_list = glob.glob(docsdir+'000*_default_norm.pha')
    for sw_path in (clock_path_list+normi_path_list):
        sw = fits.open(sw_path)
        sw_data = sw[2].data
        sw_time = (sw[2].header['TSTART']+sw[2].header['TSTOP'])/2
        fscale = (sw_time-126230400)/(12.6*365.26*86400)
        flux_sw = sw_data['FLUX']
        wave_sw = sw_data['LAMBDA']
        fig = plt.figure()
        plt.plot(wave_sw, flux_sw, 'k-',label='Swift')
        plt.plot(wave_sw, flux_sw*fscale, 'k--', alpha=0.3, label='Swift (no fscale)')
        if mod_list == None:
            pass
        else:
            for mod_single in mod_list:
                type_mod = mod_single['type']
                if type_mod == 'kurucz':
                    mod_single = normRefl(mod_single, wave_norm=3800)
                    f = interpolate.interp1d(wave_sw, flux_sw*fscale)
                    mod_single['flux'] = mod_single['flux']*f(3800)
                flux_mod = mod_single['flux']
                wave_mod = mod_single['wave']
                plt.plot(wave_mod, flux_mod,color=c_dict[type_mod],label=type_mod)
        title_list = sw_path.split('/')[-1].split('_')
        obsid = title_list[0]
        ext = title_list[1]
        if 'norm.pha' in title_list:
            mode = 'norminal'
        else:
            mode = 'clocked'
        plt.title(star+' '+obsid+' '+ext+' '+mode)
        plt.xlim(800,7000)
        plt.ylim(0,2e-13)
        plt.legend()
        plt.show()
        plt.close()

# wd1057,p177,bd25,agk81,gd153,p41,wd0320
#star_compare('p41',mod_dict)

# WD1057+719
#obsid='00055205001'
#obsid='00055211004'
#obsid='00055216001'
##obsid='00055900069'
#star='wd1057'
#obs_name = 'wd1057_719_stisnic_008.fits'
#mod_name = 'wd1057_719_mod_006.fits'

# P177
##obsid='00056760022'
#obsid='00056760024'
#obsid='00056761002'
##obsid='00056761003'
#obsid='00056762002'
#obsid='00056763002'
#obsid='00056763004'
#star = 'p177'
#obs_name = 'p177d_stisnic_008.fits'
#mod_name = None#'p177d_mod_003.fits'

# BD25
#obsid = '00057538004'
#star = 'bd25'
#obs_name = 'bd_25d4655_002.fits'
#mod_name = None

#AGK81
obsid = '00057530004'
#obsid = '00057530006'
##obsid = '00057530012'
##obsid = '00057530014'
star = 'agk81'
obs_name = 'agk_81d266_005.fits'#'agk_81d266_stisnic_007.fits'
mod_name = None


#standard_star(star, obsid, obs_name=obs_name, mod_name=mod_name)

#average_refl('juno')


'''
from scipy.ndimage import rotate
juno = fits.open('/Users/zexixing/Research/swiftASTER/data/juno/00091026003/uvot/image/sw00091026003ugu_dt.img')
img = juno[1].data
angle=35.18855991054861
if angle<40:
    theta=180-angle
else:
    theta=angle
#plt.imshow(img,vmin=0,vmax=10)
x,y=(1524.0543,585.86735)
i,j = (x-1,y-1)
new,new_point = rot(img, (i,j), theta)
plt.imshow(new,vmin=0,vmax=10)
plt.plot(new_point[0],new_point[1],'ko')
plt.show()
'''
#plot_refl()
#comparePara()
#compareAster('juno')
#compareAnalog()
#plot_refl_mean(swift=True, mmt=True, iue=True, 
#               juno=True, hebe=False, ceres=False, iris=False, massalia=False, allaster=False)
#plot_refl_sun()
#refl_smooth()

#aster_dict = {#'ceres':1, 'pallas':2, 'juno':3, 'vesta':4, 'hebe':6, 'iris':7, 'flora':8, 
              #'hygiea':10, 'eunomia':15, 'massalia':20, 'lutetia':21, 'themis':24, 'nysa':44,
              #'dembowska':349, 'eros':433, 'scheila':596, 
              #'toutatis':4179, 
#              'psyche':16,
#             }
#coord = '_1_default'
#output_name = 'motion_log2.txt'
#asterMotion(aster_dict,output_name,coord)

#swift16flux(['ceres','vesta','juno','hebe','eunomia','massalia','flora','iris','dembowska'],singleObs=False,iue=True)

#swift16flux(['iris'])
#from asjust_para import extraAster
'''
aster_list = ['ceres','juno','hebe','dembowska','eunomia',
              'flora','hygiea','iris','massalia','nysa','pallas',
              'vesta','themis','scheila','lutetia','toutatis']
for aster in aster_list:
    swift16flux([aster])
'''
#aster_list = ['toutatis','juno']
#aster_list = ['themis']
#swift16flux(aster_list,iue=True,part=False)#['00091608001','00091597001'])#'00091605001','00091600001','00091597001'])

'''
for aster in aster_list:
    idlist = ugrismlist(aster) 
    swift16flux([aster],iue=False)
'''







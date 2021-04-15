#This module determines the offset of the anchor position along the dispersion
#Inputs: aster(str), obsid(str/'mean'), sun('analog'/'colina96')
#Outputs: the offset
#1811

from tools import *
from uvotpy_support import *
from read_deal import *
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import optimize, interpolate
import numpy as np
import pandas as pd
import os

aster_list = ['toutatis','vesta','scheila','ceres',
              'massalia','nysa','pallas','lutetia',
              'juno','hebe','dembowska','eunomia',
              'flora','hygiea','iris','themis']
remove_list = [#'00091503001', #flora -- star
                '00091223001','00091225001','00091229001', #pallas -- star(compact)
                '00091198002', #'00091022002','00091197002', vesta -- star
                '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                '00091538001',#'00091540001', #nysa -- star(compact)
                '00091591001', #themis -- wrong position
                '00091207002', #juno -- star in slit
                '00091237001', #dembowska -- star in slit
                '00091268001', '00091501001', '00091503001', '00091507001', #flora -- star in slit
                '00091559001', #hygiea -- star in slit
                '00091519001', '00091521001', '00091523001', '00091525001', #iris -- star in slit
                '00091598001', '00091593001', '00091595001', #themis -- star in slit
                '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                #'00091027001', '00091241001', #dembowska -- try trackwidths
                ]
remove_list = ['00091503001', #flora -- star
                '00091223001','00091225001','00091229001', #pallas -- star(compact)
                '00091198002', #'00091022002','00091197002', vesta -- star
                '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                '00091538001','00091540001',#nysa -- star(compact)
                '00091591001', #themis -- wrong position
                '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                ]
iue_list = ['juno','ceres','iris','hebe','massalia']
mmt_list = ['juno','iris','hebe','massalia']

def check_refl(aster, sun, adjname='_smeargauss'):
    '''plot the asteroid's mean reflectance spectrum to check 
    if it is smooth enough to be fitted.

    Parameters
    ----------
    aster : str or list of str
       name of the asteroid(s), lowercase letters
    sun : str
       `analog` or `colina96`,
       solar spectrum to be used to derive reflectance
    
    kwargs : dict
    - **adjname** : str
       `_smeargauss`, `_offset`, etc
       a part of the file name which marks used paras

    Returns
    -------
    None, plot the reflectance of the inputted asteroid(s)
    Not save the figures
    '''
    # read in solar flux spectrum
    solar_spec_dict = solar_flux(sun, binbox=20)
    # loop for every inputted asteroid
    aster_spec_dict_list = {}
    if isinstance(aster,str):
        aster = [aster]
    else: pass
    aster_list = aster
    for aster in aster_list:
        # read in the asteroid's flux spectrum
        flux_dict = mean_aster_flux(aster,adjname = adjname,binbox=10)
        # derive the refl spectrum(spectra)
        refl_dict = getRefl(flux_dict, solar_spec_dict, aster=aster)
        # plot
        wave = refl_dict['wave']
        flux = refl_dict['flux']
        plt.plot(wave, flux, 'k',linewidth=0.7)
        plt.title(aster+' reflectance - quick look')
        plt.ylabel('Normalized reflectance')
        plt.xlabel('wavelength (Angstrom)')
        plt.ylim(-0.5,np.max(flux[int(len(wave)/4):])*1.5)
        plt.xlim(1600,7000)
        plt.show()

#check_refl(aster_list, 'colina96')

def fit_refl(aster, sun, adjname='_smeargauss', band=[3000,5000], method='linear', plot=False, obsid=False):
    '''fit the reflectance spectrum

    Parameters
    ----------
    aster : str or list of str
       name of the asteroid(s), lowercase letters
    sun : str
       `analog` or `colina96`,
       solar spectrum to be used to derive reflectance
    
    kwargs : dict
    - **adjname** : str
       `_smeargauss`, `_offset`, etc
       a part of the file name which marks used paras
    - **band** : list (2 elements)
       the wavelength range to be fitted
    - **method** : str
       function used to make fitting
       `linear` : y = a*x + b
       `poly` : y = a*x^2 + b*x + c
    - **plot** : bool
       whether or not to show the fitting plot

    Returns
    -------
    When *plot* = False : 
    coef : list
       Provides the fitting coefficients. 
       `linear` : [a, b]
       `poly` : [a, b, c]
    
    When *plot* = True :
       None, just plot
    '''
    # define the funtions for different methods
    if method == 'linear':
        def fun(x, a, b):
            return a*x+b
        init = [0.0005,-1.0]
    elif method == 'poly':
        def fun(x, a, b, c):
            return a*x*x+b*x+c
        init = [-0.001,3000, -3]
    # read refl
    solar_spec_dict = solar_flux(sun, binbox=20)
    if obsid:
        flux_dict = aster_flux(aster, obsid, adjname, binbox=10)
    else:
        flux_dict = mean_aster_flux(aster,adjname = adjname,binbox=10)
    refl_dict = getRefl(flux_dict, solar_spec_dict, aster=aster)  
    wave = refl_dict['wave']
    refl = refl_dict['flux']
    wave_tofit = wave[(wave>band[0])&(wave<band[1])]
    refl_tofit = refl[(wave>band[0])&(wave<band[1])]
    # fit
    coef, coef_err = optimize.curve_fit(fun,wave_tofit,refl_tofit,init)
    # plot or return
    if plot == True:
        if method == 'linear':
            a, b = coef
            plt.plot(wave, fun(wave, a, b), 'k--',linewidth=0.5)
        elif method == 'poly':
            a, b, c = coef
            plt.plot(wave, fun(x, a, b, c), 'k--',linewidth=0.5)
        plt.plot(wave, refl, 'k',linewidth=0.7)
        plt.axvline(x=band[0],c='gray',lw=1,alpha=0.2)
        plt.axvline(x=band[1],c='gray',lw=1,alpha=0.2)
        plt.ylim(-0.5,np.max(refl[int(len(wave)/4):])*1.5)
        plt.xlim(1600,7000)
        plt.ylabel('Normalized reflectance')
        plt.xlabel('wavelength (Angstrom)')
        plt.title(aster+' reflectance - fit')
        plt.show()
    else:
        return coef

#for aster in aster_list:
#    fit_refl(aster, 'colina96', plot=True)
#print(fit_refl('juno', 'colina96', plot=False))

def fit_offset(aster, obsid, sun, adjname='_smeargauss',
               band=[3000,5000],refl_method='linear',
               red_table=False, aster_dict=False, coef=False,
               output='offsetWave'):
    '''fit asteroids' flux with a reddened solar flux spectrum
    to determine offset of wavelength

    Parameters
    ----------
    aster : str or list of str
       name of the asteroid(s), lowercase letters
    obsid : str (11-digit)
       the observation to be fitted
    sun : str
       `analog` or `colina96`,
       solar spectrum to be used to derive reflectance
    
    kwargs : dict
    - **adjname** : str
       `_smeargauss`, `_offset`, etc
       a part of the file name which marks used paras
    - **band** : list (2 elements)
       the wavelength range to be fitted
    - **refl_method** : str
       function used to fit the refl's reddening
       `linear` : y = a*x + b
       `poly` : y = a*x^2 + b*x + c
    - **red_table** : dict or False
       whether or not to use a pre-determined table of reddening 
       coefficients, the inputted value will be used if yes
    - **aster_dict** : dict or bool
       whether or not to use a pre-determined asteroid's flux spectrum;
       if false, the aster_dict will be determined from the .pha data file
    - **coef** : list or bool
       whether or not to use a pre-determined offset (& norm);
       if false, the coefficients will be determined by fitting.
    - **output** : str
       `plot`, `offsetWave` or `offsetPixel`

    Returns
    -------
    When *output* = `plot` : 
    None, just plot
    
    When *output* = `offsetWave` :
    offset in unit of wavelength (A)

    When *output* = `offsetPixel` :
    offset in unit of pixel

    When *output* = `iteration` :
    '''
    # read in flux spectra
    # --- sun ---
    solar_dict = solar_flux(sun, binbox=False)
    solar_dict = normRefl(solar_dict) # normalized at 3800A by default
    wave = solar_dict['wave']
    flux = solar_dict['flux']
    wave_tofit = wave[(wave>band[0])&(wave<band[1])]
    flux_tofit = flux[(wave>band[0])&(wave<band[1])]
    # --- obs ---
    if aster_dict == False:
        aster_dict = aster_flux(aster, obsid, adjname, binbox=False) 
        aster_dict = normRefl(aster_dict) # normalized at 3800A by default
        # read or determine coefficients of the refl's reddening
        if red_table == False:
            red = fit_refl(aster, sun, adjname, band=band, method=refl_method,obsid=obsid)
        else: red = red_table
        a, b = red
    else: a, b = (0, 1)
    # define the fitting function to determine offset
    ###sun = interpolate.interp1d(solar_dict['wave'],solar_dict['flux'],fill_value="extrapolate")
    ###def fun(x, A, offset):
    ###    return A*sun(x+offset)*(a*x+b)
    aster_f = interpolate.interp1d(aster_dict['wave'],aster_dict['flux'],fill_value="extrapolate")
    def fun(x, A, offset):
        return A*aster_f(x+offset)/(a*x+b)
    init = [1,0]
    # do fitting to obtain offset
    if coef == False:
        coef, coef_err = optimize.curve_fit(fun,wave_tofit,flux_tofit,init)
        #print(coef)
    else: pass
    # plot or return
    if output == 'plot': 
        # set figure
        fig = plt.figure()
        ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4]) # left, bottom, width, height
        ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4],xticklabels=[])
        # ax1
        ax1.plot(wave,flux,'k-',linewidth=1) # sun
        ax1.plot(wave,fun(wave, coef[0], coef[1]),'r-',\
                 linewidth=0.5,label='after offset') # aster
        ax1.set_xlim(band[0],band[1])
        ax1.set_ylim(np.min(flux_tofit)-0.1,np.max(flux_tofit)*1.1)
        # ax2
        ax2.plot(wave,flux,'k-',linewidth=1) # sun
        ax2.plot(wave,fun(wave, coef[0], 0), color='b',linestyle='-',\
                 linewidth=0.5,label='before offset') # aster
        ax2.set_xlim(band[0],band[1])
        ax2.set_ylim(np.min(flux_tofit)-0.1,np.max(flux_tofit)*1.1)
        # label
        ax1.set_xlabel('wavelength (Angstrom)')
        ax1.legend()
        ax2.legend()
        if coef[1] < 0: dire = 'right'
        elif coef[1] > 0: dire = 'left'
        else: dire = 'no offset'
        plt.title(aster+' '+obsid+'  offset:'+'%.2f'%abs(coef[1])+'A ('+dire+')')
        plt.show()
    elif output == 'offsetWave':
        return coef[1]
    elif output == 'offsetPixel':
        #hdu = read_fits(aster, obsid, adjname)[2].data
        #binwidth_f = interpolate.interp1d(hdu['LAMBDA'],hdu['BINWIDTH'],fill_value="extrapolate")
        #bw_ank = binwidth_f(2600)
        #print(bw_ank)
        bw_ank = 3.2
        return coef[1]/bw_ank
    elif output == 'iteration':
        return coef, {'wave':wave,'flux':fun(wave, coef[0], coef[1])}
    else:
        pass

#for aster in aster_list:
#    obs_list = ugrismlist(aster,remove=False,filt='uvgrism')
#    for obsid in obs_list:
#        fit_offset(aster, obsid, 'colina96', band=[3000,5000],output='plot')
#fit_offset('juno', '00091026003', 'colina96', adjname='_offset', band=[3000,5000],output='plot')

def iter_offset(aster, obsid, sun, adjname='_smeargauss', output='plot'):
    '''use fit_offset to inernatively determine the final offset

    Returns
    -------
    When *output* = `plot` : 
    None, just plot
    
    When *output* = `offsetWave` :
    offset in unit of wavelength (A)

    When *output* = `offsetPixel` :
    offset in unit of pixel
    '''
    coef0, aster_dict0 = fit_offset(aster, obsid, sun, adjname=adjname, output='iteration')
    offset = coef0[1]
    norm = coef0[0]
    coef = coef0
    aster_dict = aster_dict0
    offset_list = [coef0[1]]
    loop = 1
    if abs(coef[1]<2): coef[1]=2.1 # do iteration at least for twice
    while abs(coef[1])>2 and loop < 10:
        coef, aster_dict = fit_offset(aster, obsid, sun, adjname=adjname, aster_dict=aster_dict, output='iteration')
        offset += coef[1]
        norm = norm*coef[0]
        offset_list.append(coef[1])
        loop += 1
    #print(offset,offset_list)
    if output == 'plot':
        coef, aster_dict = fit_offset(aster, obsid, sun, coef=[norm,offset], output='iteration')
        solar_dict = solar_flux(sun, binbox=False)
        solar_dict = normRefl(solar_dict) # normalized at 3800A by default
        plt.plot(aster_dict0['wave'],aster_dict0['flux'],'r-',linewidth=0.7,alpha=0.6)
        plt.plot(aster_dict['wave'],aster_dict['flux'],'k-',linewidth=0.7,alpha=0.6)
        plt.plot(solar_dict['wave'],solar_dict['flux'],'k-',linewidth=0.7)
        if offset < 0: dire = 'right'
        elif offset > 0: dire = 'left'
        plt.title(aster+' '+obsid+'  offset:'+'%.2f'%abs(offset)+'A ('+dire+')')
        plt.ylim(0,2)
        plt.xlim(3000,5000)
        plt.show()
    elif output == 'offsetWave':
        return offset
    elif output == 'offsetPixel':
        return offset/3.2
    else:
        pass

#for aster in aster_list:
#    obs_list = ugrismlist(aster,remove=False,filt='uvgrism')
#    for obsid in obs_list:
#        iter_offset(aster, obsid, 'colina96')
#iter_offset('juno', '00091026003', 'colina96', output='plot')
#iter_offset('hebe', '00091274001', 'colina96', output='plot')

def stat_offset(plot):
    '''plot statistic figurs
    
    Parameters
    ----------
    plot : str
    `hist`: plot histogram of the offset of all observations
    `motionVSoffset`: plot x-motion vs offset for all observations

    Return
    ------
    None
    '''
    offset_list = []
    offset_dict = {}
    for aster in aster_list:
        obs_list = ugrismlist(aster,remove=False,filt='uvgrism')
        for obsid in obs_list:
            offset = fit_offset(aster, obsid, 'colina96', band=[3000,5000],output='offsetWave')
            offset_list.append(offset)
            offset_dict[obsid] = offset
    # --- histogram ---
    if plot == 'hist':
        binbox = binGenerate([0,40], 5)
        plt.hist(offset_list,bins=binbox)
        plt.title('histogram of x-offset')
        plt.show()
    # --- offset vs motion
    elif plot == 'motionVSoffset':
        file_name = 'motion_log.txt'
        file_path = get_path('../docs/'+file_name)
        data = pd.read_csv(file_path, sep=' ', header=0)
        data_group = data.groupby('ASTEROID')
        motion_dict = {}
        for i in range(0,len(data)):
            motion_para = np.sqrt(data['MOTION_P'][i]**2+data['MOTION_V'][i]**2)
            obsid = '000'+str(data['OBS_ID'][i])
            aster = data['ASTEROID'][i]
            motion_dict[obsid] = (aster, motion_para)
        #offset_dict = sorted(offset_dict.items(), key=lambda d:float(d[0]))
        #motion_dict = sorted(motion_dict.items(), key=lambda d:float(d[0]))
        c_dict = {'ceres':'firebrick', 'pallas':'red', 'juno':'darksalmon', 'vesta':'gold', 'hebe':'olivedrab', 'iris':'mediumspringgreen', 'flora':'lightseagreen', 
                'hygiea':'deepskyblue', 'eunomia':'royalblue', 'massalia':'navy', 'lutetia':'blue', 'themis':'mediumpurple', 'nysa':'darkorchid',
                'dembowska':'plum', 'eros':'mediumvioletred', 'scheila':'palevioletred', 'toutatis':'sandybrown', 
                }
        for obsid in offset_dict.keys():
            label = motion_dict[obsid][0]
            plt.plot(offset_dict[obsid], motion_dict[obsid][1], 'o', \
                    color=c_dict[motion_dict[obsid][0]], MarkerSize=3)
        plt.xlabel('x-offset (A)')
        plt.ylabel('total-motion (pixel)')
        #plt.legend()
        plt.show()

#stat_offset('motionVSoffset')

def compare_offset_afbe(aster, obsid, sun, adjname_be, adjname_af):
    # sun
    '''
    solar_dict = solar_flux(sun, binbox=False)
    red = fit_refl(aster, sun, band=[3000,5000], method='linear')
    a, b = red
    solar_dict['flux'] = (a*solar_dict['wave']+b)*solar_dict['flux']
    solar_dict = normRefl(solar_dict)
    wave = solar_dict['wave']
    flux = solar_dict['flux']
    '''
    #wave_tofit = wave[(wave>band[0])&(wave<band[1])]
    #flux_tofit = flux[(wave>band[0])&(wave<band[1])]
    # aster
    flux_dict_be = aster_flux(aster, obsid, adjname_be, binbox=False)
    flux_dict_af = aster_flux(aster, obsid, adjname_af, binbox=False)
    #offset = iter_offset(aster, obsid, sun, adjname=adjname_be, output ='offsetWave')
    #flux_dict_be['wave'] = flux_dict_be['wave'] - offset
    #plt.plot(wave,flux,'k-',linewidth=0.7)
    plt.plot(flux_dict_be['wave'],flux_dict_be['flux'],'b-',linewidth=0.5)
    plt.plot(flux_dict_af['wave'],flux_dict_af['flux'],'r-',linewidth=0.5)
    plt.show()

#compare_offset_afbe('juno', '00091026003', 'colina96', '_smeargauss', '_radec')

def ank_c2anker(aster, obsid, sun, adjname='_smeargauss', output='newanker'):
    # read in y_offset, angle
    hdul = read_fits(aster, obsid, adjname)
    yoff = hdul[0].header['SPEC_OFF']
    angle = 180 - hdul[0].header['ANGLE']
    # read in anker
    anker_x = hdul[0].header['DETX_X']
    anker_y = hdul[0].header['DETY_X']
    # get x_offset
    xoff = iter_offset(aster, obsid, sun, adjname=adjname, output='offsetPixel')
    # calculate offsets of anker in det coordinate
    alpha = angle*(np.pi/180.)
    anker_xoff = xoff*np.cos(alpha) - yoff*np.sin(alpha)
    anker_yoff = xoff*np.sin(alpha) + yoff*np.cos(alpha)
    if output=='newanker':
        return anker_x+anker_xoff, anker_y+anker_yoff
    elif output=='plot':
        # plot det image
        raw_path = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'+\
                   obsid+'/uvot/image/sw'+obsid+'ugu_dt.img'
        raw = fits.open(raw_path)
        img = raw[1].data
        hdr = raw[1].header
        crpix = crpix1,crpix2 = hdr['crpix1'],hdr['crpix2'] 
        anker = np.array([anker_x, anker_y])
        anker_new = anker + np.array([anker_xoff,anker_yoff])
        ankerimg = anker - np.array([1100.5,1100.5])+crpix
        ankerimg_new = anker_new - np.array([1100.5,1100.5])+crpix 
        fig1 = plt.figure()
        plt.imshow(img, vmin=5,vmax=40,origin='lower')
        plt.plot(ankerimg[0], ankerimg[1], 'wx', MarkerSize=5)
        plt.plot(ankerimg_new[0], ankerimg_new[1], 'kx', MarkerSize=5)
        plt.show()
        # plot rotated img
        s1 = 0.5*img.shape[0]
        s2 = 0.5*img.shape[1]
        d1 = -(s1 - ankerimg[1])   # distance of anker to centre img 
        d2 = -(s2 - ankerimg[0])
        n1 = 2.*abs(d1) + img.shape[0] + 400  # extend img with 2.x the distance of anchor 
        n2 = 2.*abs(d2) + img.shape[1] + 400
        #return img, hdr, s1, s2, d1, d2, n1, n2
        if 2*int(n1/2) == int(n1): n1 = n1 + 1
        if 2*int(n2/2) == int(n2): n2 = n2 + 1
        c1 = n1 / 2 - ankerimg[1] 
        c2 = n2 / 2 - ankerimg[0]
        n1 = int(n1)
        n2 = int(n2)
        c1 = int(c1)
        c2 = int(c2)
        cval = -1.0123456789
        a = np.zeros( (n1,n2), dtype=float) + cval
        a[c1:c1+img.shape[0],c2:c2+img.shape[1]] = img
        import scipy.ndimage as ndimage
        b = ndimage.rotate(a,angle,reshape = False,order = 1,mode = 'constant',cval = cval)
        e2 = int(0.5*b.shape[0])
        c = b[e2-100:e2+100,:]
        ank_c = [ (c.shape[0]-1)/2+1, (c.shape[1]-1)/2+1 , 0, c.shape[1]]
        #from test import findBackground
        #findBackground(c,yloc_spectrum=ank_c[0])
        fig2 = plt.figure()
        plt.imshow(c,vmin=5,vmax=40,origin='lower')
        plt.plot(ank_c[1],ank_c[0],'wx',MarkerSize=5)
        plt.plot(ank_c[1]+xoff,ank_c[0]+yoff,'kx',MarkerSize=5)
        plt.show()
    else:
        pass

#print(ank_c2anker('juno', '00091026003', 'colina96','_smeargauss',output='plot'))
#ank_c2anker('dembowska', '00091239001', 'colina96','_smeargauss',output='plot')
#ank_c2anker('iris', '00091519001', 'colina96','_smeargauss',output='plot')
#ank_c2anker('pallas','00091227001','colina96','_smeargauss',output='plot')
#ank_c2anker('standard/p177', '00056763002', 'colina96','_4_default',output='plot')

def phi_fit(t,rx,ry):
   xf,yf,xp1,yp1 = getCalData(160,calfile=None,chatter=0)
   xp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), xp1 ,chatter=0)
   yp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), yp1 ,chatter=0)
   return xp1i*(1-t)/2+yp1i*(1+t)/2

def anker2sky(aster, obsid, sun, adjname, output='radec'):
    # new Xphi, Yphi
    hdul = read_fits(aster, obsid, adjname)
    Xphi = hdul[0].header['POSX_AS']
    Yphi = hdul[0].header['POSY_AS']
    anker_new = ank_c2anker(aster, obsid, sun, adjname=adjname, output='newanker')
    xyphi, xyphi_err = optimize.curve_fit(phi_fit,[-1,1],anker_new,[Xphi,Yphi])
    if output == 'xyphi':
        return xyphi
    # anker in det mm
    xyphi = np.array([xyphi[0],xyphi[1]])
    anker_uvw1det_offset = xyphi*3600/0.502
    indir='/Users/zexixing/Research/swiftASTER/data/'+aster+'/'+obsid+'/uvot/image'
    filestub = 'sw'+obsid
    specfile, lfilt1, lfilt1_ext, lfilt2, lfilt2_ext, attfile = \
        fileinfo(filestub,1,directory=indir,wheelpos=160,chatter=0)
    if type(lfilt1)==type(None):
        lfilt1 = 'fk'
    if lfilt1 == "fk" : 
        l2filter = "uvw1"
    else: l2filter = lfilt1
    anker_uvw1det = anker_uvw1det_offset + np.array(boresight(filter=l2filter,date=209952100))
    x1, y1 = (anker_uvw1det-np.array([1100.5,1100.5]))*0.009075
    # ra dec
    from astropy import wcs
    if lfilt1_ext == None: 
        lfext = 1
    else: 
        lfext = lfilt1_ext
    ffile = indir+'/'+filestub+'uw1_sk.img'
    if lfilt1 == 'wh'   : ffile = indir+'/'+filestub+'uwh_sk.img'
    if lfilt1 == 'u'    : ffile = indir+'/'+filestub+'uuu_sk.img'
    if lfilt1 == 'v'    : ffile = indir+'/'+filestub+'uvv_sk.img'
    if lfilt1 == 'b'    : ffile = indir+'/'+filestub+'ubb_sk.img'
    if lfilt1 == 'uvm2' : ffile = indir+'/'+filestub+'um2_sk.img'
    if lfilt1 == 'uvw2' : ffile = indir+'/'+filestub+'uw2_sk.img'
    if lfilt1 == 'fk'   : ffile = indir+'/'+filestub+'ufk_sk.img'
    hf = fits.getheader(ffile,lfext)
    W1 = wcs.WCS(hf,)
    W2 = wcs.WCS(hf,key='D',relax=True) 
    xpix, ypix = W2.wcs_world2pix(x1, y1, 0)
    ra, dec = W1.wcs_pix2world(xpix,ypix,0)
    ra = float(ra)
    dec = float(dec)
    if output == 'radec':
        return ra, dec
    # ra dec in hmsdms
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    c=SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
    if output == 'radec_hmsdms':
        return c.to_string('hmsdms')
    else:
        pass

#print(anker2sky('juno', '00091026003', 'colina96', '_smeargauss','radec'))
#print(anker2sky('psyche', '00014031003', 'colina96', '_1_default','radec'))

def compare_refl_afbe(aster,adjname_be,adjname_af,output='refl',wave_norm=4600):
    solar_spec_dict = solar_flux('colina96',binbox=20)
    idlist = ugrismlist(aster)
    for obsid_rm in remove_list:
        try:    idlist.remove(obsid_rm)
        except:    pass
    refl_list_be = []
    flux_list_be = []
    refl_list_af = []
    flux_list_af = []
    for obsid in idlist:
        # flux
        hdul_be = read_fits(aster, obsid, adjname_be)
        hdu_be = hdul_be[2].data
        flux_dict_be = {'wave':hdu_be['LAMBDA'],'flux':hdu_be['FLUX'],'sigma':hdu_be['FLUXERR']}
        flux_list_be.append(flux_dict_be)
        hdul_af = read_fits(aster, obsid, adjname_af)
        hdu_af = hdul_af[2].data
        flux_dict_af = {'wave':hdu_af['LAMBDA'],'flux':hdu_af['FLUX'],'sigma':hdu_af['FLUXERR']}
        flux_list_af.append(flux_dict_af)
    flux_dict_be = coadd_(flux_list_be,aster)
    flux_dict_af = coadd_(flux_list_af,aster)
    #plt.plot(flux_dict_af['wave'],flux_dict_af['flux']/flux_dict_af['sigma']) # SNR-Wave plot
    #plt.show()
    #flux_dict['wave']=flux_dict['wave']
    #return flux_dict_af
    if output == 'flux':
        # sun
        red = fit_refl(aster, 'colina96', band=[3000,5000], method='linear')
        a, b = red
        solar_spec_dict['flux'] = (a*solar_spec_dict['wave']+b)*solar_spec_dict['flux']
        f = interpolate.interp1d(flux_dict_af['wave'],flux_dict_af['flux'],fill_value="extrapolate")
        solar_spec_dict = normRefl(solar_spec_dict,wave_norm=wave_norm)
        solar_spec_dict['flux'] = solar_spec_dict['flux']*f(wave_norm)
        plt.plot(solar_spec_dict['wave'],solar_spec_dict['flux'],'k-',label='colina96',linewidth=0.5,alpha=0.4)
        plt.plot(flux_dict_be['wave'],flux_dict_be['flux'],'b-',label='before x-offset correct',linewidth=0.7)
        plt.plot(flux_dict_af['wave'],flux_dict_af['flux'],'r-',label='after x-offset correct',linewidth=0.7)
        plt.xlim(1680,6800)
        plt.ylim(0,np.max(flux_dict_af['flux'][int(len(flux_dict_af['wave'])/4):])*1.1)
        plt.fill_between([5000,6800], [0,0], [1e-11,1e-11], color='k', alpha=0.1, hatch='/')
        plt.title(aster+' - '+'averaged flux')
        plt.legend()
        plt.show()
    elif output == 'refl':
        refl_dict_be = getRefl(flux_dict_be, solar_spec_dict, aster=aster)
        refl_dict_af = getRefl(flux_dict_af, solar_spec_dict, aster=aster)
        refl_be = normRefl(refl_dict_be,wave_norm=wave_norm)
        refl_af = normRefl(refl_dict_af,wave_norm=wave_norm)
        #wave_be, flux_be = refl_be['wave'], refl_be['flux']
        #wave_af, flux_af = refl_af['wave'], refl_af['flux']
        wave_be, flux_be = refl_dict_be['wave'], refl_dict_be['flux']
        wave_af, flux_af, err_af = refl_dict_af['wave'], refl_dict_af['flux'], refl_dict_af['sigma']
        f = interpolate.interp1d(wave_af,flux_af,fill_value="extrapolate")
        plt.plot(wave_be,flux_be,'k-',label='before x-offset correct',linewidth=0.5,alpha=0.4)
        plt.plot(wave_af,flux_af,'k-',label='after x-offset correct',linewidth=0.5)
        #plt.fill_between(wave_af, flux_af-err_af, flux_af+err_af, color='k', alpha=0.2)
        plt.xlim(1680,6800)
        plt.ylim(0,np.max(flux_be[int(len(wave_be)/4):])*1.5)
        plt.fill_between([5000,6800], [0,0], [10,10], color='k', alpha=0.1, hatch='/')
        plt.title(aster+' - averaged reflectance spectrum')
        
        if aster in mmt_list:
            # mmt
            mmt_dict = {'juno':'3_Juno_27Aug09_refscl.csv',
                        'hebe':'6_Hebe_18Jun13_refscl.csv',
                        'massalia_1':'20_Massalia_UT11Dec06_refscl.csv',
                        'massalia_2':'20_Massalia_UT12Dec06_refscl.csv',
                        'iris_1':'7_Iris_11Dec06_refscl.csv',
                        'iris_2':'7_Iris_12Dec06_refscl.csv'}
            wave_norm_mmt=4600#4600
            refl_mmt = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+mmt_dict[aster], 'csv', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=',', skiprows=0)
            refl_mmt = normRefl(refl_mmt, wave_norm=wave_norm_mmt)
            plt.plot(refl_mmt['wave'],refl_mmt['flux']*(f(wave_norm_mmt)/f(wave_norm)),'r:',label='MMT',linewidth=0.5,alpha=0.5)
        if aster in iue_list:
            #iue
            wave_norm_iue=3070#3200
            refl_iue = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+aster+'_iue.txt', 'txt', 0, 1,
                                sigma_index=None, flag_index=None,
                                ext=1, sep=' ', skiprows=0)
            if aster == 'juno':
                refl_iue = normRefl({'wave':np.array(refl_iue['wave'])*10,'flux':np.array(refl_iue['flux'])}, wave_norm=wave_norm_iue)
            else:
                refl_iue = normRefl({'wave':np.array(refl_iue['wave']),'flux':np.array(refl_iue['flux'])}, wave_norm=wave_norm_iue)
            plt.plot(refl_iue['wave'],refl_iue['flux']*(f(wave_norm_iue)/f(wave_norm)),'b:',label='IUE',linewidth=0.5,alpha=0.5)
        '''
        spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass2.txt', 
                                'txt', 0, 1,  sigma_index=2,ext=1, sep='\t', skiprows=0)
        spec_dict = binSpec(np.array(spec_dict['wave'])*1e4, np.array(spec_dict['flux']), binGenerate(np.array(spec_dict['wave'])*1e4, 10))
        spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=4800)
        plt.plot(spec_dict['wave'], np.array(spec_dict['flux'])*(f(4800)/f(wave_norm)), 'g--', label='smass2',linewidth=0.5)

        spec_dict = readinSpec('/Users/zexixing/Research/swiftASTER/docs/'+aster+'/smass1.txt', 
                                'txt', 0, 1,  sigma_index=None,ext=1, sep=' ', skiprows=0)
        spec_dict = binSpec(spec_dict['wave'], np.array(spec_dict['flux']), binGenerate(spec_dict['wave'], 10), sigma=None)#np.array(spec_dict['sigma'])*norm)
        spec_dict = normRefl({'wave':np.array(spec_dict['wave']),'flux':np.array(spec_dict['flux'])},wave_norm=4800)
        plt.plot(spec_dict['wave'], np.array(spec_dict['flux'])*(f(4800)/f(wave_norm)), 'g:', label='smass1',linewidth=0.5)#spec_dict['sigma'])
        '''
        #plt.ylabel('Averaged reflectance spectra (normalized)')
        #plt.xlabel('Wavelength (A)')
        plt.legend()
        plt.show()

#for aster in ['dembowska']:#aster_list: 
#    compare_refl_afbe(aster,'_smeargauss','_radec',output='refl')

#remove_list.append('00091027001', '00091241001') #dembowska -- try trackwidths
#for aster in ['dembowska']:#aster_list: 
#    compare_refl_afbe(aster,'_radec','_tw',output='refl')
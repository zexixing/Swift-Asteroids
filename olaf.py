from tools import *
from output_plot import *
from read_deal import *
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def dt2fits(file_path, aster, obsid):
    hdul = fits.open(file_path)
    hdr = hdul[1].header[:208]
    hdr['OBJECT'] = '('+str(asterid(aster))+') '+aster.capitalize()
    hdr['PSCALE1'] = abs(float(hdul[1].header['CDELT1S'])*3600)
    hdr['PSCALE2'] = abs(float(hdul[1].header['CDELT2S'])*3600)
    #hdr['NORTH'] = 360-float(hdul[1].header['PA_PNT'])
    hdr['PRODUCTN'] = '('+str(asterid(aster))+') '+aster.capitalize()+' '+hdul[1].header['OBS_ID']+'raw ugrism detector image'
    output = fits.PrimaryHDU(hdul[1].data,header=hdr)
    outputl = fits.HDUList([output])
    output_path = get_path('../docs/olaf/'+aster+obsid+'ugu_dt'+'.fits')
    outputl.writeto(output_path)
    #os.system('mv '+output_path+' '+output_path[:-3])

#obsid = '00091026003'
#aster_name = 'juno'
#obs_name = 'sw'+obsid+'ugu_dt.img'
#obs_path = get_path('../data/'+aster_name+'/'+obsid+'/uvot/image/'+obs_name)
#dt2fits(obs_path,obsid)

def effwave():
    from read_deal import readarf, readinSpec
    f_arf, f_err = readarf(order=1, smooth=True, iferr=True, filt='uvgrism')
    solar_spec_name = '00091706004.pha'
    solar_path = '/Users/zexixing/Research/swiftASTER/docs/'+solar_spec_name
    solar_spec_dict = readinSpec(solar_path, 
                                'fits', 'LAMBDA', 'FLUX', sigma_index='FLUXERR', 
                                flag_index='QUALITY', ext=2, sep='\t', skiprows=0)
    wave = solar_spec_dict['wave']
    flux = solar_spec_dict['flux']
    return np.sum(f_arf(wave)*wave*flux)/np.sum(f_arf(wave)*flux)

#print(effwave())

def writeRefl(aster,obsid,solar_spec_dict,smooth=True,save=True,data=False):
    # flux
    adjname = '_smeargauss'
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    aster_spec_name = obsid+adjname+'.pha'
    spec_path= docsdir+aster_spec_name
    hdul = fits.open(spec_path)
    hdu = hdul[2].data
    # refl
    #solar_spec_dict = getSolarFlux(solar_spec_name, if_smooth=True)
    refl_dict = getRefl(aster_spec_name, solar_spec_dict, aster=aster)
    #output
    output_dict = {'Pixel':hdu['PIXNO'],'Wavelength':hdu['LAMBDA'],
                   'Flux':hdu['FLUX'],'Delta Flux':hdu['FLUXERR'],
                   'Reflectance':refl_dict['flux'],'Delta Reflectance':refl_dict['sigma']}

    #f = interpolate.interp1d(solar_spec_dict['wave'],solar_spec_dict['flux'],fill_value="extrapolate")
    #plt.plot(hdu['LAMBDA'],hdu['FLUX'],'r-',lw=1)
    #plt.plot(hdu['LAMBDA'],f(hdu['LAMBDA']),'b-',lw=1)
    #plt.plot(hdu['LAMBDA'],hdu['FLUX']/f(hdu['LAMBDA']),lw=1)
    #plt.show()
    if data:
        return output_dict
    output_df = pd.DataFrame.from_dict(output_dict)
    columns = [('Pixel',''),('Wavelength','[Angstrom]'),
               ('Flux','[erg/cm2/s/A]'),('Delta Flux','[erg/cm2/s/A]'),
               ('Reflectance','[]'),('Delta Reflectance','[]')]
    output_df.columns=pd.MultiIndex.from_tuples(columns)
    if save:
        output_df.to_csv('/Users/zexixing/Research/swiftASTER/docs/olaf/'+aster+'_'+obsid+'.csv',
                         index=False,header=False)#float_format='%.4f')
    else:
        data = output_df.to_string(index=False,header=False)
        data = data[1:]
        data = data.replace('\n   ','\n')
        data = data.replace('\n  ','\n')
        data = data.replace('\n ','\n')
        data = data.replace('     ',',')
        data = data.replace('    ',',')
        data = data.replace('   ',',')
        data = data.replace('  ',',')
        data = data.replace(' ',',')
        data = data.replace(',,',',')
        return data
    #plt.plot(hdu['LAMBDA'],refl_dict['flux'])
    #plt.show()

def averageRefl(aster,solar_spec_dict,save=True,data=False):
    # get ids
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
        # refl
        aster_spec_name = obsid+adjname+'.pha'
        refl_dict = getRefl(aster_spec_name, solar_spec_dict, aster=aster)
        refl_list.append(refl_dict)
        # flux
        spec_path= docsdir+aster_spec_name
        hdul = fits.open(spec_path)
        hdu = hdul[2].data
        flux_dict = {'wave':hdu['LAMBDA'],'flux':hdu['FLUX'],'sigma':hdu['FLUXERR']}
        flux_list.append(flux_dict)
    refl_dict = coadd_(refl_list,aster)
    flux_dict = coadd_(flux_list,aster)
    #output
    output_dict = {'Wavelength':refl_dict['wave'],
                   'Flux':flux_dict['flux'],'Delta Flux':flux_dict['sigma'],
                   'Reflectance':refl_dict['flux'],'Delta Reflectance':refl_dict['sigma']}
    if data:
        return output_dict
    output_df = pd.DataFrame.from_dict(output_dict)
    columns = [('Wavelength','[Angstrom]'),
               ('Flux','[erg/cm2/s/A]'),('Delta Flux','[erg/cm2/s/A]'),
               ('Reflectance','[]'),('Delta Reflectance','[]')]
    output_df.columns=pd.MultiIndex.from_tuples(columns)
    if save:
        output_df.to_csv('/Users/zexixing/Research/swiftASTER/docs/olaf/'+aster+'_average'+'.csv',
                         index=False,columns=False)
    else:
        data = output_df.to_string(index=False,header=False)
        data = data[1:]
        data = data.replace('\n   ','\n')
        data = data.replace('\n  ','\n')
        data = data.replace('\n ','\n')
        data = data.replace('     ',',')
        data = data.replace('    ',',')
        data = data.replace('   ',',')
        data = data.replace('  ',',')
        data = data.replace(' ',',')
        data = data.replace(',,',',')
        return data
    #plt.plot(refl_dict['wave'],flux_dict['flux'])
    #plt.show()

def obslog(aster,save=True):
    from astropy.time import Time
    output_name = aster+'_obs-log.csv'
    output_path = '/Users/zexixing/Research/swiftASTER/docs/olaf/'+output_name
    datadir = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'
    idlist = ugrismlist(aster,filt='allfilter')
    obstable=''
    tstart_list = []
    tstop_list = []
    for obsid in idlist:
        data_path = datadir+obsid+'/uvot/image/'
        try:
            obs_name = 'sw'+obsid+'ugu_dt.img' ####
            obs_path = data_path+obs_name
            hdul = fits.open(obs_path)
            hdr = hdul[1].header
            exp = hdr['EXPOSURE']
        except:
            filelist = os.listdir(data_path)
            filelist = [i for i in filelist if 'sk.img' in i]
            obs_name = filelist[0]
            obs_path = data_path+obs_name
            hdul = fits.open(obs_path)
            hdr = hdul[0].header
            ext_num = len(hdul)
            exp = 0
            for i in range(1,ext_num):
                exp += float(hdul[i].header['EXPOSURE'])
        tstart = hdr['DATE-OBS']
        tstop = hdr['DATE-END']
        tstart_list.append(tstart)
        tstop_list.append(tstop)
        filt = hdr['FILTER']
        tmid = (Time(tstop)-Time(tstart))/2+Time(tstart)
        obj = Horizons(id=asterid(aster),location='@swift',epochs=tmid.jd)
        eph = obj.ephemerides()[0]
        rh = eph['r']
        rhv = eph['r_rate']
        delta = eph['delta']
        ra = eph['RA']
        dec = eph['DEC']
        phase = eph['alpha']
        c = SkyCoord(ra, dec, frame='icrs', unit='deg')
        c = c.to_string('hmsdms')
        ra, dec = c.split(' ')
        obsinfo = obsid+','+tstart+','+tstop+','+tmid.isot+','\
                  +str(exp)+','+filt+','+str(rh)+','+str(rhv)+','\
                  +str(delta)+','+str(phase)+','+ra+','+dec+'\n'
        obstable += obsinfo
    if save:
        title = 'OBS_ID,START_TIME,END_TIME,MIDTIME,EXP_TIME,\
                FILTER,HELIO.,HELIOV.,DELTA,PHASE,RA,DEC\n'
        units = ',[UT],[UT],[UT],[s],,[AU],[km/s],[AU],[deg],[hms],[dms]\n'
        f = open(output_path, 'w')
        f.write(title+units+obstable)
        f.close()
    else:
        return obstable, tstart_list[0],tstop_list[-1]

def sunFlux():
    from scipy import interpolate
    aster='solar_analog'
    obsid = '00091706004'
    #obs_name = 'sw'+obsid+'ugu_dt.img'
    #file_path = get_path('../data/'+aster+'/'+obsid+'/uvot/image/'+obs_name)
    #hdul = fits.open(file_path)
    pha_path = get_path('../docs/'+obsid+'.pha')
    hdul = fits.open(pha_path)
    pha = hdul[2].data
    smooth = coadd_([{'wave':pha['LAMBDA'],'flux':pha['FLUX'],'sigma':pha['FLUXERR']}],
                     aster='solar_analog')
    f = interpolate.interp1d(smooth['wave'],smooth['flux'],fill_value="extrapolate")
    e = interpolate.interp1d(smooth['wave'],smooth['sigma'],fill_value="extrapolate")
    output_dict = {'Pixel':pha['PIXNO'],'Wavelength':pha['LAMBDA'],
                   'Flux':pha['FLUX'],'Delta Flux':pha['FLUXERR'],
                   'Smoothed Flux':f(pha['LAMBDA']),'Smoothed Delta Flux':e(pha['LAMBDA'])}
    output_df = pd.DataFrame.from_dict(output_dict)
    data = output_df.to_string(index=False,header=False)
    data = data[1:]
    data = data.replace('\n   ','\n')
    data = data.replace('\n  ','\n')
    data = data.replace('\n ','\n')
    data = data.replace('     ',',')
    data = data.replace('    ',',')
    data = data.replace('   ',',')
    data = data.replace('  ',',')
    data = data.replace(' ',',')
    data = data.replace(',,',',')
    return data, hdul[1].header['DATE-OBS'], hdul[1].header['DATE-END']

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
                '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                #'00091027001', '00091241001', #dembowska -- try trackwidths
                ]
                
def olafCsv(aster, obsid, optype):
    if optype=='individual':
        obs_name = 'sw'+obsid+'ugu_dt.img'
        file_path = get_path('../data/'+aster+'/'+obsid+'/uvot/image/'+obs_name)
        hdul = fits.open(file_path)
        data = writeRefl(aster,obsid,solar_spec_dict,smooth=True,save=False)
        div = ',,,,'
        level = 'Calibrated'
        columns = 'Pixel,Wavelength,Flux,Delta Flux,Reflectance,Delta Reflectance'
        datatype = 'integer,real,real,real,real,real'
        units = ',Angstrom,erg/cm2/s/A,erg/cm2/s/A,,'
        missing = '-999,-999,-999,-999,-999,-999'
        output_name = aster+'_'+obsid+'.csv'
        product_name = '('+str(asterid(aster))+') '+aster.capitalize()+' '+obsid+' individual spectrum'
        start_time = hdul[0].header['DATE-OBS']
        stop_time = hdul[0].header['DATE-END']
        target_name = '('+str(asterid(aster))+') '+aster.capitalize()
        descrip = 'Flux spectrum and normalized reflectance spectrum (normalized at 3800A) for an individual exposure of '+\
                  str(asterid(aster))+aster.capitalize()+'; observation ID: '+obsid+'.'
        data_descrip = 'Pixel number with origin at anchor,\
                        1-order wavelength corresponding to the pixel number,\
                        Calibrated 1-order flux,\
                        1-sigma errors of the flux,\
                        Reflectance from the calibrated flux,\
                        1-sigma error of the reflectance'
    elif optype=='average':
        data = averageRefl(aster,solar_spec_dict,save=False)
        div = ',,,'
        level = 'Derived'
        columns = 'Wavelength,Flux,Delta Flux,Reflectance,Delta Reflectance'
        datatype = 'real,real,real,real,real'
        units = 'Angstrom,erg/cm2/s/A,erg/cm2/s/A,,'
        missing = '-999,-999,-999,-999,-999'
        output_name = aster+'_average.csv'
        product_name = '('+str(asterid(aster))+') '+aster.capitalize()+' average spectrum'
        idlist = [int(obsid) for obsid in ugrismlist(aster,remove=remove_list)]
        idstr = ' '.join(ugrismlist(aster,remove=remove_list))
        first_obs = min(idlist)
        last_obs = max(idlist)
        obs_name = 'sw'+'000'+str(first_obs)+'ugu_dt.img'
        file_path = get_path('../data/'+aster+'/'+'000'+str(first_obs)+'/uvot/image/'+obs_name)
        hdul = fits.open(file_path)
        start_time = hdul[0].header['DATE-OBS']
        obs_name = 'sw'+'000'+str(last_obs)+'ugu_dt.img'
        file_path = get_path('../data/'+aster+'/'+'000'+str(last_obs)+'/uvot/image/'+obs_name)
        hdul = fits.open(file_path)
        stop_time = hdul[0].header['DATE-END']
        target_name = '('+str(asterid(aster))+') '+aster.capitalize()
        descrip = 'Averaged flux spectrum and averaged normalized reflectance spectrum (normalized at 3800A) for '+\
                  str(asterid(aster))+aster.capitalize()+\
                  '; ID of averaged observations: '+idstr+'.'
        data_descrip = '1-order wavelength corresponding to the pixel number,\
                        Calibrated 1-order flux,\
                        1-sigma errors of the flux,\
                        Reflectance from the calibrated flux,\
                        1-sigma error of the reflectance'
    elif optype == 'log':
        div = ',,,,,,,,,,'
        level = 'Partially Processed'
        data, start_time, stop_time= obslog(aster,save=False)
        columns = 'OBS_ID,START_TIME,END_TIME,MIDTIME,EXP_TIME,\
                FILTER,HELIO.,HELIOV.,DELTA,PHASE,RA,DEC'
        datatype = 'string,datetime,datetime,datetime,real,\
                   string,real,real,real,real,string,string'
        units = ',UT,UT,UT,s,,AU,km/s,AU,deg,hms,dms'
        missing = '"Missing","0000-00-00T00:00:00","0000-00-00T00:00:00","0000-00-00T00:00:00",\
                  -999,"Missing",-999,-999,-999,-999,"Missing","Missing"'
        output_name = aster+'_obs-log.csv'
        product_name = '('+str(asterid(aster))+') '+aster.capitalize()+' observing log'
        target_name = '('+str(asterid(aster))+') '+aster.capitalize()
        descrip = 'Observing inforation and orbital information (from JPL/Horizons) for this observation campaign of'+\
                  str(asterid(aster))+aster.capitalize()
        data_descrip = 'Swift observing ID,\
                        Start time of the observation,\
                        End time of the observation,\
                        Midtime of the observation,\
                        Total exposure time of the observation,\
                        Filter used in the observation,\
                        Heliocentric distance of the object at midtime (JPL/Horizons),\
                        Heliocentric velocity of the object at midtime (JPL/Horizons),\
                        Apparent range of the object relative to Swift at midtime (JPL/Horizons),\
                        The Sun-Target-Observer angle at midtime (JPL/Horizons),\
                        RA of the object relative to Swift at midtime (JPL/Horizons),\
                        DEC of the object relative to Swift at midtime (JPL/Horizons)'
    elif optype == 'sun': 
        data, start_time, stop_time = sunFlux()
        div = ',,,,'
        level = 'Calibrated'
        columns = 'Pixel,Wavelength,Flux,Delta Flux,Smoothed Flux,Smoothed Delta Flux'
        datatype = 'integer,real,real,real,real,real'
        units = ',Angstrom,erg/cm2/s/A,erg/cm2/s/A,erg/cm2/s/A,erg/cm2/s/A'
        missing = '-999,-999,-999,-999,-999,-999'
        output_name = 'solar_analog.csv'
        product_name = 'solar analog - BD124536 spectrum'
        target_name = 'BD-124536'
        descrip = 'Flux spectrum of BD-124536, a G2V-type solar analog, observed by Swift uvgrism and calibrated by UVOTPY. The spectrum is binned by 10A to smooth, and interpolated to derive asteroids\' reflectance spectra.'
        data_descrip = 'Pixel number with origin at anchor,\
                        1-order wavelength corresponding to the pixel number,\
                        Calibrated 1-order flux,\
                        1-sigma errors of the flux,\
                        Smoothed calibrated 1-order flux,\
                        1-sigma errors of the smoothed flux'
    output_path = '/Users/zexixing/Research/swiftASTER/docs/olaf/'+output_name
    f = open(output_path,'w')
    f.write('# Keywords\n'\
            + 'File Name,'+output_name+div+'\n'\
            + 'Product Name,'+product_name+div+'\n'\
            + 'Product Description,"'+descrip+'"'+div+'\n'\
            + 'Start Time,'+start_time+div+'\n'\
            + 'Stop Time,'+stop_time+div+'\n'\
            + 'Target Name,'+target_name+div+'\n'\
            + 'Target Type,Asteroid'+div+'\n'\
            + 'Author List,"Bodewits, D.; Xing, Z."'+div+'\n'\
            + 'Product Processing Level,'+level+div+'\n'\
            + 'Science Search Facet,"Spectral Image,Grayscale,Lightcurve,Tabulated,Physical Properties"'+div+'\n'\
            + 'Product Wavelength Ranges,Ultraviolet'+div+'\n'\
            + 'Observing System Bookmark,"uvot"'+div+'\n'\
            + '# Column Definitions'+div+'\n'\
            + columns+'\n'\
            + datatype+'\n'\
            + units+'\n'\
            + missing+'\n'\
            + data_descrip+'\n'\
            + '# Data\n'\
            + data)
    f.close()

def graphsObs(aster,obsid,solar_spec_dict):
    data = writeRefl(aster,obsid,solar_spec_dict,smooth=True,save=False,data=True)
    flux = data['Flux']
    flux_err = data['Delta Flux']
    refl = data['Reflectance']
    refl_err = data['Delta Reflectance']
    wave = data['Wavelength']
    fig = plt.figure()
    ax21 = plt.subplot(2,1,1)
    plt.plot(wave,flux,color='darkred',lw=1)
    plt.fill_between(wave, flux-flux_err, flux+flux_err, color='darkred', alpha=0.2)
    qf = np.max(flux[int(len(wave)*0.3):int(len(wave)*0.7)])
    plt.ylim(0.001*qf,1.2*qf)
    plt.xlim(1663,6800)
    plt.ylabel('1st order flux\n('+r'$erg~cm^{-2}s^{-1}\AA^{-1}$'+')')
    ax22 = plt.subplot(2,1,2)
    plt.plot(wave,refl,color='darkred',lw=1)
    plt.fill_between(wave, refl-refl_err, refl+refl_err, color='darkred', alpha=0.2)
    qf = np.max(refl[int(len(wave)*0.3):int(len(wave)*0.7)])
    plt.ylim(0.001*qf,1.2*qf)
    plt.xlim(1663,6800)
    plt.xlabel('Wavelength('+r'$\AA$'+')')
    plt.ylabel('Reflectance\n(Normalized at 3800'+r'$\AA$'+')')
    output_name = 'refl_'+obsid+'.png'
    output_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+output_name
    plt.savefig(output_path,dpi=300)
    plt.close()

def graphsAve(aster,solar_spec_dict):
    data = averageRefl(aster,solar_spec_dict,save=False,data=True)
    flux = np.array(data['Flux'])
    flux_err = np.array(data['Delta Flux'])
    refl = np.array(data['Reflectance'])
    refl_err = np.array(data['Delta Reflectance'])
    wave = np.array(data['Wavelength'])
    # flux
    fig = plt.figure()
    plt.plot(wave,flux,color='darkred',lw=1)
    plt.fill_between(wave, flux-flux_err, flux+flux_err, color='darkred', alpha=0.2)
    qf = np.max(flux[int(len(wave)*0.3):int(len(wave)*0.7)])
    plt.ylim(0.001*qf,1.2*qf)
    plt.xlim(1663,6800)
    plt.ylabel('1st order averaged flux ('+r'$erg~cm^{-2}s^{-1}\AA^{-1}$'+')')
    output_name = 'average_flux_'+aster+'.png'
    output_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+output_name
    plt.savefig(output_path,dpi=300)
    plt.close()
    # refl
    fig = plt.figure()
    plt.plot(wave,refl,color='darkred',lw=1)
    plt.fill_between(wave, refl-refl_err, refl+refl_err, color='darkred', alpha=0.2)
    qf = np.max(refl[int(len(wave)*0.3):int(len(wave)*0.7)])
    plt.ylim(0.001*qf,1.2*qf)
    plt.xlim(1663,6800)
    plt.xlabel('Wavelength('+r'$\AA$'+')')
    plt.ylabel('Averaged reflectance (Normalized at 3800'+r'$\AA$'+')')
    output_name = 'average_refl_'+aster+'.png'
    output_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'+output_name
    plt.savefig(output_path,dpi=300)
    plt.close()


solar_spec_dict = getSolarFlux('00091706004.pha', if_smooth=True)
aster = 'juno'
idlist = ugrismlist(aster)
#graphsAve(aster,solar_spec_dict)
for obsid in idlist:
    print(aster+' '+obsid)
    obs_name = 'sw'+obsid+'ugu_dt.img'
    obs_path = get_path('../data/'+aster+'/'+obsid+'/uvot/image/'+obs_name)
    #graphsObs(aster,obsid,solar_spec_dict)
    dt2fits(obs_path, aster, obsid)
    #olafCsv(aster, obsid, 'individual')
#olafCsv(aster, None, 'average')


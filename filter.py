from aper_phot import *
from tools import *
from read_deal import readarf, ugrismlist, coadd_
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.time import Time
from _mypath import thisdir
from scipy import interpolate
from scipy import optimize
from ank_x_offset import fit_refl

aster_list = ['vesta','scheila','ceres',
              'massalia','nysa','pallas','lutetia',
              'juno','hebe','dembowska','eunomia',
              'flora','hygiea','iris','themis']
remove_list = [ '00091223001','00091225001','00091229001', #pallas -- star(compact)
                #'00091198002', #'00091022002','00091197002', vesta -- star
                '00091220001','00091216001', #'00091218001', lutetia -- star(compact)
                '00091540001','00091538001', #nysa -- star(compact)
                '00091591001', #themis -- wrong position
                '00091207002', #juno -- star in slit
                '00091237001', #dembowska -- star in slit
                '00091268001', '00091503001', '00091507001', #flora -- star in slit, @'00091501001'
                #@'00091559001', #hygiea -- star in slit
                #@'00091519001', '00091521001', '00091523001', '00091525001', #iris -- star in slit
                '00091598001', #@'00091593001', '00091595001', #themis -- star in slit
                '00091600001', '00091609001', '00091606001', '00091597001', '00091603001', #toutatis -- no filter
                #@'00091027001', #dembowska -- star in slit
                '00091242001', #ceres -- star in slit
                '00091022002', #vesta -- star
                ]
rmstar_list = ['00091525001','00091519001','00091523001', '00091521001',#iris
               '00091532001', #nysa
               '00091593001','00091595001',#themis
               '00091027001',#dembowska
               '00091505001','00091501001',#flora
               '00091559001',#hygiea
               ]
c_list = ['firebrick', 'red', 'darksalmon', 'gold', 'olivedrab', 'mediumspringgreen', 'lightseagreen', 
          'deepskyblue', 'royalblue', 'navy', 'blue', 'mediumpurple', 'darkorchid',
          'plum', 'mediumvioletred', 'palevioletred', 'sandybrown']
def generate_c():
    for c in c_list:
        yield c
color_g=generate_c()

def coi_corr(raw):
    a1 = 0.066
    a2 = -0.091
    a3 = 0.029
    a4 = 0.031
    ft = 0.0110329
    x = raw*ft
    y = 1+a1*x+a2*np.power(x,2)+\
        a3*np.power(x,3)+a4*np.power(x,4)
    alp = 0.9842
    z = -np.log(1-alp*raw*ft)/(alp*ft)
    f = y*z#/raw
    #f = np.nan_to_num(f,nan=1.0,posinf=1.0,neginf=1.0)
    return f

def radec_horizons(obs_log_name,horizons_id,filt):
    obs_log_path = get_path('../docs/'+obs_log_name)
    f = open(obs_log_path,'r')
    skiprows=1
    obs_dic = {}
    for line in f.readlines()[skiprows:]:
        line = line.strip()
        line = line.split(',')
        obs_id = line[0]
        filt_key = line[5]
        if filt_key == filt:
            img_name = 'sw000'+obs_id+filt_dict[filt_key]+'_sk.img'+file_dict[filt]
            relative_path = '../data/'+obj_name+'/000'+obs_id+'/uvot/image/'
            img_path = get_path('../docs/'+relative_path+img_name)
            hdul = fits.open(img_path)
            ext_header = hdul[1].header
            w = WCS(ext_header)
            midtime = float(line[3])
            obj = Horizons(id=horizons_id, location='@swift', epochs=float(line[3]))
            eph = obj.ephemerides()[0]
            ra = eph["RA"]
            dec = eph["DEC"]
            px, py = w.wcs_world2pix(ra, dec, 1)
            obs_dic[obs_id] = {'midtime':midtime,
                                'ra':ra,'dec':dec,
                                'px':px,'py':py}
    return obs_dic

def lightcurve(obj_name,filt):
    obs_dic = radec_horizons(obs_log_name,horizons_id,filt)
    mag_list = []
    mag_err_list = []
    jdtime_list = []
    for obs_id in obs_dic.keys():
        src_cen = [obs_dic[obs_id]['px'],obs_dic[obs_id]['py']]
        img_name = 'sw000'+obs_id+filt_dict[filt]+'_sk.img'+file_dict[filt]
        relative_path = '../data/'+obj_name+'/000'+obs_id+'/uvot/image/'
        result = aper_phot(img_name, filt_dict[filt],
                           src_cen, src_r=30,
                           bg_center=src_cen, bg_r = [250,300],
                           src_method='total_mean', bg_method='donut_mean',
                           step=5, mask_img = False, 
                           start = False, relative_path=relative_path,ext=1)
        mag_list.append(result[2][0])
        mag_err_list.append(result[2][1])
        jdtime_list.append(obs_dic[obs_id]['midtime'])
    plt.plot(jdtime_list, mag_list, marker='o',label=filt)
    plt.title(obj_name)
    plt.xlabel('JD time (days)')
    plt.ylabel('mag')
    plt.legend()
'''
filt_dict = {'UVW2':'uw2','UVM2':'um2'}
file_dict = {'UVW2':'','UVM2':'.gz'}
obs_log_name = 'obs-log_hebe.csv'
horizons_id = 6
obj_name = 'hebe'
lightcurve(obj_name,'UVW2')
lightcurve(obj_name,'UVM2')
plt.show()
'''
def check(aster, spec):
    from read_deal import binSpec, binGenerate
    idlist = ugrismlist(aster)
    adjname = '_rmstar'
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    flux_list = []
    arf = readarf(filt='um2')
    for obsid in idlist:
        if spec == 'bkg':
            if obsid in rmstar_list:
                aster_spec_name = obsid+'_rmstar_back.pha'
            else:
                aster_spec_name = obsid+'_radec_back.pha'
            spec_path= docsdir+aster_spec_name
            hdul = fits.open(spec_path)
            hdu = hdul[1].data
            plt.plot(hdu['CHANNEL'][::-1],hdu['BKG_COUNTS']/hdu['EXPOSURE'],label=obsid)
        elif spec == 'flux':
            if obsid in rmstar_list:
                aster_spec_name = obsid+'_rmstar.pha'
            else:
                aster_spec_name = obsid+'_radec.pha'
            spec_path= docsdir+aster_spec_name
            hdul = fits.open(spec_path)
            data = hdul[2].data
            flux = data['FLUX']
            flux[:200]=0
            wave=data['LAMBDA']
            #data = binSpec(wave, flux, binGenerate([1200,7000],30), sigma=None, flag=None)
            #flux = data['flux']
            #wave = data['wave']
            hnu = 6.626e-27*2.9979e10*1e8/wave
            ea_filt = arf(wave)
            #plt.plot(wave,flux*ea_filt/hnu,label=obsid)
            plt.plot(wave,flux,label=obsid)
            plt.xlim(1800,3000)
    plt.title(aster)
    plt.legend()
    plt.show()


#check('juno','flux')
#check('juno','bkg')
    
def average_spec(aster):
    idlist = ugrismlist(aster,remove='rmstar')
    adjname = '_rmstar'
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    flux_list = []
    for obsid in idlist:
        # flux
        if obsid in rmstar_list:
            aster_spec_name = obsid+'_rmstar'+'.pha'
        else:
            aster_spec_name = obsid+'_radec'+'.pha'
        print(aster_spec_name)
        spec_path= docsdir+aster_spec_name
        hdul = fits.open(spec_path)
        hdu = hdul[2].data
        flux_dict = {'wave':hdu['LAMBDA'],'flux':hdu['FLUX'],'sigma':hdu['FLUXERR']}
        flux_list.append(flux_dict)
    flux_dict = coadd_(flux_list,aster,binbox=10)
    plt.plot(flux_dict['wave'],flux_dict['flux'])
    plt.show()
    return flux_dict
#average_spec('dembowska')

#red_a, red_b = fit_refl('iris', 'colina96', adjname='_rmstar', band=[2900,3900])
#print(red_a,red_b)
def read_sun():
    from read_deal import getSolarFlux
    data = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    #sun_flux = data['flux']*(red_a*data['wave']+red_b)
    f = interpolate.interp1d(data['wave'],data['flux'])
    #f = interpolate.interp1d(data['wave'],sun_flux)
    return f
sun_f = read_sun()

def sun_replac_fitting(wave, a, b, c):
    return a*sun_f(wave)*(b*wave+c)

def sun_replace(wave, flux, sigma):
    # cut out ： 2600 - 3000 for the asteroid & sun
    cut_dw = 2500
    cut_up = 2850
    wave_2fit = wave[wave>cut_dw]
    flux_2fit = flux[wave>cut_dw]
    sigma_2fit = sigma[wave>cut_dw]
    flux_2fit = flux_2fit[wave_2fit<cut_up]
    sigma_2fit = sigma_2fit[wave_2fit<cut_up]
    wave_2fit = wave_2fit[wave_2fit<cut_up]
    # do fitting & obtain the a
    coef, coef_err = optimize.curve_fit(sun_replac_fitting,wave_2fit,flux_2fit,[1e-15, 0.0005, -0.9])
    a,b,c=tuple(coef)
    a_err,b_err,c_err = tuple(np.sqrt(np.diag(coef_err)))
    # replace <2650
    wave_rpl = 2650
    wave_2rpl = wave[wave<wave_rpl]
    flux_2rpl = a*sun_f(wave_2rpl)*(b*wave_2rpl+c)
    #sigma_2rpl = a_err*sun_f(wave_2rpl)*np.sqrt(np.power(b_err*wave_2rpl,2)+np.power(c_err,2))
    sigma_2rpl = flux_2rpl*2
    flux_stay = flux[wave>=wave_rpl]
    sigma_stay = sigma[wave>wave_rpl]
    flux_new = np.append(flux_2rpl,flux_stay)
    sigma_new = np.append(sigma_2rpl,sigma_stay)
    sigma_new = sigma
    #plt.plot(wave,flux_new)
    #plt.errorbar(wave,flux_new,yerr=sigma_new)
    #plt.show()
    return flux_new,sigma_new,a,b,c


def predcount(aster,asterid,sp_name,len_id,filt):
    # sp_name -> sp_path, read sp, read arf
    rela_path = '/Users/zexixing/Research/swiftASTER/docs/'+aster+'/'
    sp_path = rela_path + sp_name
    data = fits.open(sp_path)[2].data
    hdr = fits.open(sp_path)[2].header
    swifttime = hdr['TSTART']
    pix = data['PIXNO']
    wave = data['LAMBDA']
    netrate = data['NETRATE']
    bkgrate = data['BGRATE1']
    flux = data['FLUX']
    sigma = data['FLUXERR']
    fscale =(swifttime-126230400.000) / (12.6*365.26*86400)
    flux = flux*fscale
    #flux[:200]=0
    flux,sigma,a,b,c = sun_replace(wave, flux, sigma)
    #if filt == 'um2':
    #    color=next(color_g)
    #    plt.plot(wave,data['FLUX']*fscale,alpha=0.3,c=color)
    #    plt.plot(wave,a*sun_f(wave)*(b*wave+c),c=color,label=sp_name[:11])
    #    #lt.plot(wave,flux,c=color,label=sp_name[:11])
    #    plt.xlim(1400,3200)
    #    plt.ylim(-0.5e-13,2e-13)
    #    plt.legend()
    #    plt.show()
    apercorr = data['APERCOR1']
    exposure = data['EXPOSURE']
    fcoi = data['SP1_COIF']
    bgcoi = data['BG1_COIF']
    senscorr = data['SENSCORR']
    binwidth = data['BINWIDTH']
    ea = data['EFFAREA1']
    arf = readarf(filt=filt)
    ea_filt = arf(wave)
    # predict
    hnu = 6.626e-27*2.9979e10*1e8/wave
    sprate_filt = flux*ea_filt*binwidth/hnu
    sprate_filt_err = sigma*ea_filt*binwidth/hnu
    #plt.plot(wave,flux*ea_filt/hnu)
    #plt.title(sp_name[:11])
    #plt.ylim(0,0.012)
    #plt.xlim(1700,7000)
    #plt.show()
    #####sprate_obs = sprate_filt#*(2/3)#/(senscorr)
    pred = np.sum(sprate_filt)
    pred_err = np.sqrt(np.sum(np.power(sprate_filt_err,2)))
    # actual
    try:
        img_name = 'sw'+len_id+filt+'_sk.img.gz'
        relative_path = '../data/'+aster+'/'+len_id+'/uvot/image/'
        hdul = fits.open(thisdir+'/'+relative_path+img_name)
    except:
        img_name = 'sw'+len_id+filt+'_sk.img'
        relative_path = '../data/'+aster+'/'+len_id+'/uvot/image/'
        hdul = fits.open(thisdir+'/'+relative_path+img_name) 
    ext_header = hdul[1].header
    w = WCS(ext_header)
    start = ext_header['DATE-OBS']
    end = ext_header['DATE-END']
    midtime = (Time(end)-Time(start))/2+Time(start)
    obj = Horizons(id=asterid, location='@swift', epochs=midtime.jd)
    eph = obj.ephemerides()[0]
    ra = eph["RA"]
    dec = eph["DEC"]
    px, py = w.wcs_world2pix(ra, dec, 1)
    print(ra,dec,px,py)
    src_cen=[px,py]
    result = aper_phot(img_name, filt,
                        src_cen, src_r=10,
                        bg_center=src_cen, bg_r = [30,70],
                        src_method='total_mean', bg_method='donut_mean',
                        step=5, mask_img = False, 
                        start = False, relative_path=relative_path,ext=1)
    actu = result[0][0]
    actu_err = result[0][1]
    print(actu,pred)
    actu = coi_corr(actu)
    pred = 1.05*pred
    print(actu,pred)
    return pred,actu,pred_err,actu_err

#pred,actu=predcount('juno',3,'00091026003_radec.pha','00091026003','uw2')
#pred,actu=predcount('hebe',6,'00091274001_2order_4width.pha','00091274001','uw2')
#pred,actu=predcount('hygiea',10,'00091550001_radec.pha','00091550001','uw2')
#pred,actu=predcount('dembowska',349,'00091027001_2order_4width.pha','00091027001','uw2')
#pred,actu=predcount('flora',8,'00091505001_2order_4width.pha','00091505001','uw2')
#pred,actu=predcount('iris',3,'00091519001_rmstar.pha','00091519001','um2') # neg
#pred,actu=predcount('iris',3,'00091521001_rmstar.pha','00091521001','um2') # pos
#pred,actu=predcount('iris',3,'00091523001_rmstar.pha','00091523001','um2') # pos/neg
#pred,actu=predcount('iris',3,'00091525001_rmstar.pha','00091525001','um2') #pos
#print(pred,actu)


aster_dict = {'juno':3, 'vesta':4, 'hebe':6, 'iris':7, 'flora':8, 
              'hygiea':10, 'eunomia':15, 'themis':24,
              #'nysa':44, 'lutetia':21, 'massalia':20, 'pallas':2,
              'dembowska':349, #'eros':433, 'scheila':596, 
              'toutatis':4179, 
              #'ceres':1,
             }
c_dict = {'ceres':'firebrick', 'pallas':'red', 'juno':'darksalmon', 'vesta':'gold', 'hebe':'olivedrab', 'iris':'mediumspringgreen', 'flora':'lightseagreen', 
          'hygiea':'deepskyblue', 'eunomia':'royalblue', 'massalia':'navy', 'lutetia':'blue', 'themis':'mediumpurple', 'nysa':'darkorchid',
          'dembowska':'plum', 'eros':'mediumvioletred', 'scheila':'palevioletred', 'toutatis':'sandybrown', 
          'wd0320':'firebrick', 'wd1057':'red', 'wd1657':'darksalmon', 'gd153':'gold', 
          'p177':'olivedrab', 'p41':'mediumspringgreen', 'agk81':'lightseagreen', 'bd25':'deepskyblue',
         }
sky_dict = {'wd0320':(50.56166667, -53.75458333),
            'wd1057':(165.1425, 71.63441667),
            'wd1657':(254.71291667, 34.31486111),
            'gd153':(194.25958333, 22.03130556),
            'p177':(239.80666667, 47.61163889),
            'p41':(222.99166667, 71.7215),
            'agk81':(140.33, 81.72433333),
            #'bd25':(237.99958333, 32.94841667),
            }
sp_dict = {'wd0320':'wd0320_539_mod_001.fits',#'wd0320_539_stis_005.fits',
           'wd1057':'wd1057_719_mod_006.fits',
           'wd1657':'wd1657_343_mod_007.fits',
           'gd153':'gd153_mod_011.fits',
           'p177':'p177d_stisnic_008.fits',#'p177d_mod_003.fits',
           'p41':'p041c_mod_003.fits',#'p041c_stisnic_008.fits',
           'agk81':'agk_81d266_stisnic_007.fits',
           'bd25':'bd_25d4655_002.fits',}

um2_dict = {'p177':['00056762001','00056763003'],
            'bd25':['00057538003']}
def fun(x, a):
    return a*x

def pre_obs(aster_dict, adjname):
    filt2 = {'l':'uw2','u':'UVW2'}
    #filt2 = {'l':'uw1','u':'UVW1'}
    aster_dict_uvm2 = {}
    aster_dict_uvw2 = {}
    for aster in aster_dict.keys():
        print(aster)
        aster_pred_list_uvm2 = []
        aster_pred_err_list_uvm2 = []
        aster_meas_list_uvm2 = []
        aster_meas_err_list_uvm2 = []
        aster_pred_list_uvw2 = []
        aster_pred_err_list_uvw2 = []
        aster_meas_list_uvw2 = []
        aster_meas_err_list_uvw2 = []
        uvg_list = ugrismlist(aster,remove=False, filt='uvgrism')
        um2_list = ugrismlist(aster,remove=False, filt='um2')
        uw2_list = ugrismlist(aster,remove=False, filt=filt2['l'])
        if aster == 'toutatis':
            for obsid_rm in remove_list:
                try:    uvg_list.remove(obsid_rm)
                except:    pass
            um2_list = uvg_list
        if len(uvg_list) == len(um2_list):
            for i in range(len(uvg_list)):
                if uvg_list[i] in remove_list:
                    break
                print(uvg_list[i])
                sp_name = uvg_list[i]+adjname+'.pha'
                if uvg_list[i] in rmstar_list:
                    sp_name = uvg_list[i]+'_rmstar.pha'
                cnt_pred, cnt_meas, cnt_pred_err, cnt_meas_err = predcount(aster,aster_dict[aster],sp_name,um2_list[i],'um2')
                aster_pred_list_uvm2.append(cnt_pred)
                aster_pred_err_list_uvm2.append(cnt_pred_err)
                aster_meas_list_uvm2.append(cnt_meas)
                aster_meas_err_list_uvm2.append(cnt_meas_err)
                print(cnt_pred, cnt_meas)
            #plt.xlim(1400,3200)
            #plt.ylim(-0.5e-13,2e-13)
            #plt.legend()
            #plt.show()
        if len(aster_pred_list_uvm2) != 0:
            aster_dict_uvm2[aster] = [aster_pred_list_uvm2,aster_meas_list_uvm2,aster_pred_err_list_uvm2,aster_meas_err_list_uvm2]
        if len(uvg_list) == len(uw2_list):
            for i in range(len(uvg_list)):
                if uvg_list[i] in remove_list:
                    break
                print(uvg_list[i])
                sp_name = uvg_list[i]+adjname+'.pha'
                if uvg_list[i] in rmstar_list:
                    sp_name = uvg_list[i]+'_rmstar.pha'
                cnt_pred, cnt_meas, cnt_pred_err, cnt_meas_err = predcount(aster,aster_dict[aster],sp_name,uw2_list[i],filt2['l'])
                aster_pred_list_uvw2.append(cnt_pred)
                aster_pred_err_list_uvw2.append(cnt_pred_err)
                aster_meas_list_uvw2.append(cnt_meas)
                aster_meas_err_list_uvw2.append(cnt_meas_err)
        if len(aster_pred_list_uvw2) != 0:
            aster_dict_uvw2[aster] = [aster_pred_list_uvw2,aster_meas_list_uvw2,aster_pred_err_list_uvw2,aster_meas_err_list_uvw2]
    fig_uvm2 = plt.figure()
    plt.title('UVM2')
    plt.xlabel('grism (counts/s)')
    plt.ylabel('filter (counts/s)')
    #plt.xlim(xmin=0)
    #plt.ylim(ymin=0)
    uvm2_pred = []
    uvm2_meas = []
    for aster in aster_dict_uvm2.keys():
        plt.scatter(aster_dict_uvm2[aster][0],aster_dict_uvm2[aster][1],color=c_dict[aster],
                    marker='*',label=aster)
        plt.errorbar(aster_dict_uvm2[aster][0],aster_dict_uvm2[aster][1],color=c_dict[aster],
                     xerr=aster_dict_uvm2[aster][2],yerr=aster_dict_uvm2[aster][3],fmt='none')
        uvm2_pred.extend(aster_dict_uvm2[aster][0])
        uvm2_meas.extend(aster_dict_uvm2[aster][1])
    if uvm2_pred:
        coef_uvm2, coef_uvm2_err = optimize.curve_fit(fun,uvm2_pred,uvm2_meas,[1])
        print('coef uvm2: ',coef_uvm2[0], coef_uvm2_err[0])
        plt.plot([0,80],[0,coef_uvm2[0]*80],'k-',lw=1)
    plt.legend()
    fig_uvw2 = plt.figure()
    plt.title(filt2['u'])
    plt.xlabel('grism (counts/s)')
    plt.ylabel('filter (counts/s)')
    #plt.xlim(xmin=0)
    #plt.ylim(ymin=0)
    uvw2_pred = []
    uvw2_meas = []
    for aster in aster_dict_uvw2.keys():
        plt.scatter(aster_dict_uvw2[aster][0],aster_dict_uvw2[aster][1],
                    color=c_dict[aster],marker='o',label=aster)
        #plt.errorbar(aster_dict_uvw2[aster][0],aster_dict_uvw2[aster][1],color=c_dict[aster],
        #             xerr=aster_dict_uvw2[aster][2],yerr=aster_dict_uvw2[aster][3],fmt='none')
        uvw2_pred.extend(aster_dict_uvw2[aster][0])
        uvw2_meas.extend(aster_dict_uvw2[aster][1])
    if uvw2_pred:
        print(uvw2_pred,uvw2_meas)
        coef_uvw2, coef_uvw2_err = optimize.curve_fit(fun,uvw2_pred,uvw2_meas,[1])
        print('coef '+filt2['u']+': ',coef_uvw2[0], coef_uvw2_err[0])
        plt.plot([0,220],[0,coef_uvw2[0]*220],'k-',lw=1)
    plt.legend()
    plt.show()

#pre_obs(aster_dict, '_radec')
#pre_obs({'juno':3}, '_radec')

def pred_standard(star,filt):
    # -----pred-----
    from read_deal import binSpec, binGenerate
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
    name_sp = sp_dict[star]
    path_sp = docsdir + name_sp
    sp = fits.open(path_sp)
    data_sp = sp[1].data
    flux_sp = data_sp['FLUX']
    wave_sp = data_sp['WAVELENGTH']
    #from read_deal import getSolarFlux
    #data_sp = getSolarFlux('sun_ref_colina96.asc.txt', if_smooth=True)
    #wave_sp = data_sp['wave'][100:]
    #flux_sp = data_sp['flux'][100:]
    data_sp = binSpec(wave_sp, flux_sp, binGenerate([1200,7000],10), sigma=None, flag=None)
    flux_sp = data_sp['flux']
    wave_sp = data_sp['wave']
    binwidth=10
    #binwidth = (wave_sp[2:]+wave_sp[1:-1])/2 - (wave_sp[1:-1]+wave_sp[:-2])/2
    #wave_sp = wave_sp[1:-1]
    #flux_sp = flux_sp[1:-1]
    arf = readarf(filt=filt)
    ea_filt = arf(wave_sp)
    hnu = 6.626e-27*2.9979e10*1e8/wave_sp
    sprate_filt = flux_sp*ea_filt*binwidth/hnu
    #plt.plot(wave_sp,flux_sp)
    #plt.show()
    pred = np.sum(sprate_filt)
    #print(pred)
    return pred
#pred_standard('p41','uw2')
def actu_standard(star,len_id,filt):
    #-----actu-----
    relative_path = '../data/standard/'+star+'/'+len_id+'/uvot/image/'
    try:
        img_name = 'sw'+len_id+filt+'_sk.img.gz'
        hdul = fits.open(thisdir+'/'+relative_path+img_name)
    except:
        img_name = 'sw'+len_id+filt+'_sk.img'
        hdul = fits.open(thisdir+'/'+relative_path+img_name) 
    #print(img_name)
    ext_header = hdul[1].header
    swifttime = (ext_header['TSTART']+ext_header['TSTOP'])/2
    printtime = ext_header['DATE-OBS']
    print(printtime)
    w = WCS(ext_header)
    ra, dec = sky_dict[star]
    px, py = w.wcs_world2pix(ra, dec, 1)
    print(ra,dec,px,py)
    src_cen=[px,py]
    result = aper_phot(img_name, filt,
                        src_cen, src_r=10,
                        bg_center=src_cen, bg_r = [30,70],
                        src_method='total_mean', bg_method='donut_mean',
                        step=5, mask_img = False, 
                        start = False, relative_path=relative_path,ext=1)
    actu = result[0][0]
    return actu

def pre_obs_standard(filt):
    import glob
    star_list = list(sky_dict.keys())
    if filt == 'uw1':
        star_list.remove('gd153')
        #star_list = ['p177']
        pass
    elif filt == 'um2':
        star_list = ['p177']
    all_pred = []
    all_actu = []
    #star_filt_dict = {}
    for star in star_list:
        print(star)
        docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+star+'/'
        datadir = '/Users/zexixing/Research/swiftASTER/data/standard/'+star+'/'
        clock_path_list = glob.glob(docsdir+'000*_default.pha')
        if star=='p177' and filt=='um2':
            clock_path_list = [docsdir+'00056762001_1_default.pha',docsdir+'00056763003_1_default.pha']
        else:
            pass
        name_sp = sp_dict[star]
        pred = pred_standard(star,filt)
        star_actu_filt = []
        star_pred_filt = []
        for clock_path in clock_path_list:
            file_name = clock_path.split('/')[-1]
            if file_name.split('_')[1] == '1':
                obsid = file_name.split('_')[0]
                print(obsid)
                try:
                    actu = actu_standard(star,obsid,filt)
                    actu=coi_corr(actu)
                    star_actu_filt.append(actu)
                    star_pred_filt.append(pred)
                    print(actu,pred)
                except:
                    pass
            else:
                pass
        all_pred.extend(star_pred_filt)
        all_actu.extend(star_actu_filt)
        #star_filt_dict[star] = [star_actu_filt, star_pred_filt]
        if star_actu_filt:
            plt.scatter(star_pred_filt,star_actu_filt,color=c_dict[star],marker='*',label=star)
    from scipy import optimize
    coef, coef_err = optimize.curve_fit(fun,all_pred,all_actu,[1])
    print('coef '+filt+': ',coef[0], coef_err[0])
    plt.plot([0,200],[0,coef[0]*200],'k-',lw=1)
    plt.ylabel('filter (counts/s)')
    plt.xlabel('predict (counts/s)')
    plt.title(filt)
    plt.legend()
    plt.show()
        
#pre_obs_standard('uw1')



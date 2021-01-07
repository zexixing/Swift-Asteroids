from aper_phot import *
from tools import *
from read_deal import readarf, ugrismlist
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.time import Time
from _mypath import thisdir

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
    flux[:200]=0
    apercorr = data['APERCOR1']
    exposure = data['EXPOSURE']
    fcoi = data['SP1_COIF']
    bgcoi = data['BG1_COIF']
    senscorr = data['SENSCORR']
    binwidth = data['BINWIDTH']
    ea = data['EFFAREA1']
    arf = readarf(filt=filt)
    ea_filt = arf(wave)
    fscale =(swifttime-126230400.000) / (12.6*365.26*86400)
    flux = flux#*fscale
    # predict
    hnu = 6.626e-27*2.9979e10*1e8/wave
    sprate_filt = flux*ea_filt*binwidth/hnu
    sprate_obs = sprate_filt#*(2/3)#/(senscorr)
    pred = np.sum(sprate_obs)
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
    src_cen=[px,py]

    result = aper_phot(img_name, filt,
                        src_cen, src_r=30,
                        bg_center=src_cen, bg_r = [250,300],
                        src_method='total_mean', bg_method='donut_mean',
                        step=5, mask_img = False, 
                        start = False, relative_path=relative_path,ext=1)
    actu = result[0][0]
    return pred,actu

#pred,actu=predcount('juno',3,'00091026003_2order_4width.pha','00091026003','uw2')
#pred,actu=predcount('hebe',6,'00091274001_2order_4width.pha','00091274001','uw2')
#pred,actu=predcount('hygiea',10,'00091550001_radec.pha','00091550001','uw2')
#pred,actu=predcount('dembowska',349,'00091027001_2order_4width.pha','00091027001','uw2')
#pred,actu=predcount('flora',8,'00091505001_2order_4width.pha','00091505001','uw2')

#print(pred,actu)

aster_list = ['vesta','scheila','ceres',
              'massalia','nysa','pallas','lutetia',
              'juno','hebe','dembowska','eunomia',
              'flora','hygiea','iris','themis']
aster_dict = {'juno':3, 'vesta':4, 'hebe':6, 'iris':7, 'flora':8, 
              'hygiea':10, 'eunomia':15, 'themis':24,
              #'nysa':44, 'lutetia':21, 'massalia':20, 'pallas':2,
              'dembowska':349, #'eros':433, 'scheila':596, 
              'toutatis':4179, 
              #'ceres':1,
             }
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
               '00091559001',#'hygiea'
               ]
c_dict = {'ceres':'firebrick', 'pallas':'red', 'juno':'darksalmon', 'vesta':'gold', 'hebe':'olivedrab', 'iris':'mediumspringgreen', 'flora':'lightseagreen', 
          'hygiea':'deepskyblue', 'eunomia':'royalblue', 'massalia':'navy', 'lutetia':'blue', 'themis':'mediumpurple', 'nysa':'darkorchid',
          'dembowska':'plum', 'eros':'mediumvioletred', 'scheila':'palevioletred', 'toutatis':'sandybrown', 
         }

def fun(x, a):
    return a*x

def pre_obs(aster_dict, adjname):
    from scipy import optimize
    aster_dict_uvm2 = {}
    aster_dict_uvw2 = {}
    for aster in aster_dict.keys():
        print(aster)
        aster_pred_list_uvm2 = []
        aster_meas_list_uvm2 = []
        aster_pred_list_uvw2 = []
        aster_meas_list_uvw2 = []
        uvg_list = ugrismlist(aster,remove=False, filt='uvgrism')
        um2_list = ugrismlist(aster,remove=False, filt='um2')
        uw2_list = ugrismlist(aster,remove=False, filt='uw2')
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
                cnt_pred, cnt_meas = predcount(aster,aster_dict[aster],sp_name,um2_list[i],'um2')
                aster_pred_list_uvm2.append(cnt_pred)
                aster_meas_list_uvm2.append(cnt_meas)
        if len(aster_pred_list_uvm2) != 0:
            aster_dict_uvm2[aster] = [aster_pred_list_uvm2,aster_meas_list_uvm2]
        if len(uvg_list) == len(uw2_list):
            for i in range(len(uvg_list)):
                if uvg_list[i] in remove_list:
                    break
                print(uvg_list[i])
                sp_name = uvg_list[i]+adjname+'.pha'
                if uvg_list[i] in rmstar_list:
                    sp_name = uvg_list[i]+'_rmstar.pha'
                cnt_pred, cnt_meas = predcount(aster,aster_dict[aster],sp_name,uw2_list[i],'uw2')
                aster_pred_list_uvw2.append(cnt_pred)
                aster_meas_list_uvw2.append(cnt_meas)
        if len(aster_pred_list_uvw2) != 0:
            aster_dict_uvw2[aster] = [aster_pred_list_uvw2,aster_meas_list_uvw2]
    fig_uvm2 = plt.figure()
    plt.title('UVM2')
    plt.xlabel('grism (counts/s)')
    plt.ylabel('filter (counts/s)')
    #plt.xlim(xmin=0)
    #plt.ylim(ymin=0)
    uvm2_pred = []
    uvm2_meas = []
    for aster in aster_dict_uvm2.keys():
        plt.scatter(aster_dict_uvm2[aster][0],aster_dict_uvm2[aster][1],color=c_dict[aster],marker='*',label=aster)
        uvm2_pred.extend(aster_dict_uvm2[aster][0])
        uvm2_meas.extend(aster_dict_uvm2[aster][1])
    if uvm2_pred:
        coef_uvm2, coef_uvm2_err = optimize.curve_fit(fun,uvm2_pred,uvm2_meas,[1])
        print('coef uvm2: ',coef_uvm2[0], coef_uvm2_err[0])
        plt.plot([0,80],[0,coef_uvm2[0]*80],'k-',lw=1)
    plt.legend()
    fig_uvw2 = plt.figure()
    plt.title('UVW2')
    plt.xlabel('grism (counts/s)')
    plt.ylabel('filter (counts/s)')
    #plt.xlim(xmin=0)
    #plt.ylim(ymin=0)
    uvw2_pred = []
    uvw2_meas = []
    for aster in aster_dict_uvw2.keys():
        plt.scatter(aster_dict_uvw2[aster][0],aster_dict_uvw2[aster][1],color=c_dict[aster],marker='o',label=aster)
        uvw2_pred.extend(aster_dict_uvw2[aster][0])
        uvw2_meas.extend(aster_dict_uvw2[aster][1])
    if uvw2_pred:
        print(uvw2_pred,uvw2_meas)
        coef_uvw2, coef_uvw2_err = optimize.curve_fit(fun,uvw2_pred,uvw2_meas,[1])
        print('coef uvw2: ',coef_uvw2[0], coef_uvw2_err[0])
        plt.plot([0,220],[0,coef_uvw2[0]*220],'k-',lw=1)
    plt.legend()
    plt.show()

pre_obs(aster_dict, '_radec')
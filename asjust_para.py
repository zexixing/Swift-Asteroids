from uvotpy import uvotgetspec
from tools import getRadec
from output_plot import ugrismlist
from ank_x_offset import fit_offset, iter_offset, anker2sky
import os

def adjustParas():
    adjname = '_radec2' # '_xxx'
    docsdir = '/Users/zexixing/Research/swiftASTER/docs'
    datadir = '/Users/zexixing/Research/swiftASTER/data/solar_analog/00091706004/uvot/image'
    os.chdir(docsdir)
    ra, dec = 248.125833, -12.763611
    #ra,dec = 248.1275598,-12.7624883
    ext = 1
    obsid = '00091706004'
    uvotgetspec.getSpec(ra,dec,obsid,ext,indir=datadir,wr_outfile=True,predict2nd=False,plot_img=True,plot_raw=True,plot_spec=True,
                        fit_second=False,chatter=1,wheelpos=160,clobber=True,highlight=False,skip_field_src=True
                        )#background_lower=[70,90],background_upper=[50,70])
    os.system('rm *.log *.ub1 *.spec')
    os.system('mv '+datadir+'/'+obsid+'_map.png'+' '+docsdir+'/'+obsid+'_map'+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_count.png'+' '+docsdir+'/'+obsid+'_count'+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_flux.png'+' '+docsdir+'/'+obsid+'_flux'+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_fit.png'+' '+docsdir+'/'+obsid+'_fit'+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_wing.png'+' '+docsdir+'/'+obsid+'_wing'+adjname+'.png')
    os.system('mv '+docsdir+'/sw'+obsid+'ugu_1ord_1_f.pha'+' '+docsdir+'/'+obsid+adjname+'.pha')

#adjustParas()

def extraAster_point(obs,aster,ra,dec,mode,ext=1):
    from tools import getRadec
    from astropy.io import fits
    from subtract import replace_img
    adjname = '_'+str(int(ext))+'_default'#'_smeargauss'#_4width' # '_xxx'
    if mode == '200':
        adjname = adjname+'_norm'
    obsid = obs# 00091204002 91205002 91206002 91207002
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/standard/'+aster
    datadir = '/Users/zexixing/Research/swiftASTER/data/standard/'+aster+'/'+obsid+'/uvot/image'
    os.chdir(docsdir)
    #ra, dec = getRadec(obsid,horizon_id,'@swift',aster)
    ##file_name = 'sw'+obsid+'ugu_dt.img'
    ##ra = fits.open(datadir+'/'+file_name)[0].header['RA_OBJ']
    ##dec = fits.open(datadir+'/'+file_name)[0].header['DEC_OBJ']
    #file_path='/Users/zexixing/Research/swiftASTER/docs/motion_log.txt'
    #x_offset = iter_offset(aster, obsid, 'colina96', output='offsetPixel')
    #x_offset = fit_offset(aster, obsid, 'colina96', adjname='_smeargauss',output='offsetPixel')
    replace = replace_img(obsid)
    uvotgetspec.getSpec(ra,dec,obsid,ext,indir=datadir,wr_outfile=True,predict2nd=True,plot_img=True,plot_raw=True,plot_spec=True,
                        fit_second=True,chatter=2,wheelpos=160,clobber=True,highlight=False,skip_field_src=True, 
                        ifmotion=False,lfilt1_ext=1,lfilt2_ext=1)
                        #background_lower=[50,75],background_upper=[50,75])
    os.system('rm *.log *.ub1 *.spec')
    os.system('mv '+datadir+'/'+obsid+'_map.png'+' '+docsdir+'/'+'map_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_count.png'+' '+docsdir+'/'+'count_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_flux.png'+' '+docsdir+'/'+'flux_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_fit.png'+' '+docsdir+'/'+'fit_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_wing.png'+' '+docsdir+'/'+'wing_'+obsid+adjname+'.png')
    os.system('mv '+docsdir+'/sw'+obsid+'ugu_1ord_'+str(int(ext))+'_f.pha'+' '+docsdir+'/'+obsid+adjname+'.pha')
    os.system('mv '+docsdir+'/sw'+obsid+'ugu_1ord_'+str(int(ext))+'_f_back.pha'+' '+docsdir+'/'+obsid+adjname+'_back.pha')

sky_dict = {'wd1057':(165.1425,71.63441667),
            'p177':(239.80666667,47.61163889),
            'bd25':(329.92491667,26.43261111),
            'agk81':(140.33,81.72433333),
            'gd153':(194.25958333,22.03130556),
            'p41':(222.99166667, 71.7215),
            'wd0320':(50.56166667, -53.75458333),
            'wd1657':(254.71291667, 34.31486111),
            'wr1':(10.86833333, 64.75983333),
            'wr4':(40.29875, 56.73047222),
            'wr52':(199.61666667, -58.13711111),
            'wr86':(259.59625, -34.4085),
            'wr121':(281.055, -3.79938889)}
obs_dict = {#'wd1057':['00055216001','00055205001','00055211004','00055900069'],
            #'p177':['00056760022','00056760024','00056761002','00056761003','00056762002','00056763002','00056763004'],
            #'bd25':['00057538004'],
            #'agk81':['00057530004','00057530006','00057530012','00057530014'],
            #'gd153':['00055505002','00055505003'],
            #'p41':['00057955002','00057964002','00057967002'],
            #'wd0320':['00054250023','00054250027','00054250029','00054255003','00054257002'],
            #'wr1':['00037905003'],
            #'wr4':['00056900002','00056900004'],
            #'wr52':['00057927003','00057928004'],
            #'wr86':['00057012002','00057015002'],
            #'wr121':['00057500003'],
            'wd1657':['00055900069']}

from astropy.io import fits
from os.path import isfile
for star in obs_dict.keys():
    print(star)
    ra,dec = sky_dict[star]
    obs_list = obs_dict[star]
    for obsid in obs_list:
        datadir = '/Users/zexixing/Research/swiftASTER/data/standard/'+star+'/'+obsid+'/uvot/image/'
        filename = 'sw'+obsid+'ugu_dt.img'
        print(obsid)
        if not isfile(datadir+filename):
            filename = 'sw'+obsid+'ugu_dt.img.gz'
            if not isfile(datadir+filename):
                continue
        file_path = datadir+filename
        mode = fits.getheader(file_path,0)['DETNAM']
        if (mode != '160') & (mode != '200'):
            continue
        exts = len(fits.open(file_path))-1
        for ext in range(1,exts+1):
            extraAster_point(obsid,star,ra,dec,mode,ext)

#DETNAM

def extraAster(obs,aster,horizon_id):
    from tools import getRadec
    from astropy.io import fits
    from subtract import replace_img
    adjname = '_rmstar'#'_smeargauss'#_4width' # '_xxx'
    obsid = obs# 00091204002 91205002 91206002 91207002
    print(aster,obsid)
    docsdir = '/Users/zexixing/Research/swiftASTER/docs/'+aster
    datadir = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'+obsid+'/uvot/image'
    os.chdir(docsdir)
    #ra, dec = getRadec(obsid,horizon_id,'@swift',aster)
    ##file_name = 'sw'+obsid+'ugu_dt.img'
    ##ra = fits.open(datadir+'/'+file_name)[0].header['RA_OBJ']
    ##dec = fits.open(datadir+'/'+file_name)[0].header['DEC_OBJ']
    ra,dec=anker2sky(aster, obsid, 'colina96', '_smeargauss', output='radec')
    print(ra,dec)
    ext = 1
    file_path='/Users/zexixing/Research/swiftASTER/docs/motion_log.txt'
    #x_offset = iter_offset(aster, obsid, 'colina96', output='offsetPixel')
    #x_offset = fit_offset(aster, obsid, 'colina96', adjname='_smeargauss',output='offsetPixel')
    replace = replace_img(obsid)
    uvotgetspec.getSpec(ra,dec,obsid,ext,indir=datadir,wr_outfile=True,predict2nd=True,plot_img=True,plot_raw=True,plot_spec=True,
                        fit_second=True,chatter=2,wheelpos=160,clobber=True,highlight=False,skip_field_src=True, 
                        ifmotion=True, motion_file=file_path,ank_c_0offset=True,replace=replace)
                        #background_lower=[50,75],background_upper=[50,75])
    os.system('rm *.log *.ub1 *.spec')
    os.system('mv '+datadir+'/'+obsid+'_map.png'+' '+docsdir+'/'+'map1_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_count.png'+' '+docsdir+'/'+'count1_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_flux.png'+' '+docsdir+'/'+'flux1_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_fit.png'+' '+docsdir+'/'+'fit1_'+obsid+adjname+'.png')
    os.system('mv '+datadir+'/'+obsid+'_wing.png'+' '+docsdir+'/'+'wing1_'+obsid+adjname+'.png')
    os.system('mv '+docsdir+'/sw'+obsid+'ugu_1ord_1_f.pha'+' '+docsdir+'/'+obsid+adjname+'.pha')
    os.system('mv '+docsdir+'/sw'+obsid+'ugu_1ord_1_f_back.pha'+' '+docsdir+'/'+obsid+adjname+'_back.pha')
    #os.system('rm *.png')
'''
# juno
aster = 'juno'
for obs in ['00091026003','00091203002','00091204002','00091205002','00091206002','00091207002']:
    extraAster(obs,aster,3)

# hebe
aster = 'hebe'
print(aster)
for obs in ['00091274001','00091510001','00091512001','00091514001','00091516001']:
    extraAster(obs,aster,6)

# ceres
aster = 'ceres'
print(aster)
for obs in ['00091242001','00091244001','00091246001','00091248001']:
    extraAster(obs,aster,1)

# iris
aster = 'iris'
print(aster)
for obs in ['00091519001','00091521001','00091523001','00091525001','00091527001']:
    extraAster(obs,aster,7)

# massalia
aster = 'massalia'
print(aster)
for obs in ['00091541001','00091543001','00091545001','00091547001','00091549001']:
    extraAster(obs,aster,20)
'''

aster_dict = {'ceres':1, 'juno':3, 'hebe':6, 'dembowska':349,
              'eunomia':15, 'flora':8, 'hygiea':10, #'eros':433, 
              'iris':7, 'lutetia':21, 'massalia':20, 'nysa':44,
              'pallas':2, 'themis':24, 'vesta':4,
              'scheila':596, 'yu':55,
              'toutatis':4179}
#aster_dict = {'juno':3}
#aster = 'juno'
#aster_id = 3
#aster_dict = {'solar_analog':1}
'''
for aster in aster_dict.keys():
    aster_id = aster_dict[aster]
    idlist = ugrismlist(aster)
    print(idlist)
    for obsid in idlist:
        extraAster(obsid,aster,aster_id)
'''

#extraAster('00091207002','juno',3)
#extraAster('00091501001','flora',8)
#extraAster('00091242001','ceres',1)
#extraAster('00091519001','iris',7)
#extraAster('00091239001','dembowska',349)
#extraAster('00091234001','eunomia',15)
#extraAster('00091609001','toutatis',4179)
#extraAster('00091706004','solar_analog',None)

#extraAster('00091525001','iris',7)
#extraAster('00091519001','iris',7)
#extraAster('00091523001','iris',7)
#extraAster('00091532001','nysa',44)
#extraAster('00091593001','themis',24)
#extraAster('00091595001','themis',24)
#extraAster('00091027001','dembowska',349)
#extraAster('00091505001','flora',8)
#extraAster('00091559001','hygiea',10)
#extraAster('00091521001','iris',7)

#extraAster('00091503001','flora',8)
#extraAster('00091507001','flora',8)
#extraAster('00091268001','flora',8)
#extraAster('00091207002','juno',3)
#extraAster('00091501001','flora',8)
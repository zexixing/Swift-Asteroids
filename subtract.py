from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from read_deal import ugrismlist
from ank_x_offset import anker2sky
from astropy import wcs
from tools import *
import pandas as pd

def test():
    rela_path = '/Users/zexixing/Research/swiftASTER/data/hygiea/'
    obsid1 = '00091554001'
    hdu1 = fits.open(rela_path+obsid1+'/uvot/image/sw'+obsid1+'ugu_dt.img')
    data1 = hdu1[1].data
    exp1 = hdu1[1].header['EXPOSURE']
    data1 = data1/exp1
    obsid2 = '00091556001'
    hdu2 = fits.open(rela_path+obsid2+'/uvot/image/sw'+obsid2+'ugu_dt.img')
    data2 = hdu2[1].data
    exp2 = hdu2[1].header['EXPOSURE']
    data2 = data2/exp2
    obsid3 = '00091559001'
    hdu3 = fits.open(rela_path+obsid3+'/uvot/image/sw'+obsid3+'ugu_dt.img')
    data3 = hdu3[1].data
    exp3 = hdu3[1].header['EXPOSURE']
    data3 = data3/exp3
    obsid4 = '00091550001'
    hdu4 = fits.open(rela_path+obsid4+'/uvot/image/sw'+obsid4+'ugu_dt.img')
    data4 = hdu4[1].data
    exp4 = hdu4[1].header['EXPOSURE']
    data4 = data4/exp4
    obsid5 = '00091552001'
    hdu5 = fits.open(rela_path+obsid5+'/uvot/image/sw'+obsid5+'ugu_dt.img')
    data5 = hdu5[1].data
    exp5 = hdu5[1].header['EXPOSURE']
    data5 = data5/exp5
    med = np.median([data1,data2,data3,data4,data5],axis=0)
    sub4 = data4-med
    plt.imshow(sub4,vmin=0,vmax=0.1)
    plt.show()

def read_data(aster, obsid):
    docsdir = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'+obsid+'/uvot/image/'
    try:
        aster_spec_name = 'sw'+obsid+'ugu_dt.img'
        spec_path= docsdir+aster_spec_name
        hdul = fits.open(spec_path)
    except:
        aster_spec_name = 'sw'+obsid+'ugu_dt.img.gz'
        spec_path= docsdir+aster_spec_name
        hdul = fits.open(spec_path)
    return hdul

def bkg_med(aster, remove_list=False):
    obs_list = ugrismlist(aster,remove=remove_list,filt='uvgrism')
    med = []
    for obsid in obs_list:
        hdu = read_data(aster, obsid)
        bkg_cr = hdu[1].data/hdu[1].header['EXPOSURE']
        med.append(bkg_cr)
    med = np.median(med,axis=0)
    return med


def shift_star(aster, main_id, star_id, offset=(0,0), fill=False, fill_sed=False):
    # ra1, dec1 & pixel1 of the aster in the star image
    #ra,dec=anker2sky(aster, star_id, 'colina96', '_smeargauss', output='radec')
    hdu_star = read_data(aster, star_id)
    exp_star = hdu_star[1].header['EXPOSURE']
    ra = hdu_star[1].header['RA_OBJ']
    dec = hdu_star[1].header['DEC_OBJ']
    data_star = hdu_star[1].data
    wS_star =wcs.WCS(header=hdu_star[1].header,key='S',relax=True,)
    ypix_star, xpix_star = wS_star.wcs_world2pix(ra,dec,0)
    # determine pixel3 of ra1, dec1 in the main image
    hdu_main = read_data(aster, main_id)
    exp_main = hdu_main[1].header['EXPOSURE']
    data_main = hdu_main[1].data
    wS_main =wcs.WCS(header=hdu_main[1].header,key='S',relax=True,)
    ypix_main, xpix_main = wS_main.wcs_world2pix(ra,dec,0)
    # median bkg
    med = bkg_med(aster)
    star_rmbkg =  data_star/exp_star - med
    star_rmbkg[np.where(star_rmbkg<0.005)] = 0
    #star_rmbkg[np.where(star_rmbkg>=0.005)] = 9999
    # shift pixel1 to pixel3 
    shift = np.zeros(data_star.shape)
    x, y = data_star.shape
    rg_offset, up_offset = offset #(-1,2)
    rg = round(abs(xpix_main-xpix_star))+rg_offset
    up = round(abs(ypix_main-ypix_star))+up_offset
    if ((xpix_main-xpix_star)>=0) & ((ypix_main-ypix_star)>=0):
        shift[rg:,up:] = star_rmbkg[:x-rg,:y-up]
    elif ((xpix_main-xpix_star)<0) & ((ypix_main-ypix_star)<0):
        shift[:x-rg,:y-up] = star_rmbkg[rg:,up:]
    elif ((xpix_main-xpix_star)<0) & ((ypix_main-ypix_star)>=0):
        print('3')
        shift[:x-rg,up:] = star_rmbkg[rg:,:y-up]
    elif ((xpix_main-xpix_star)>=0) & ((ypix_main-ypix_star)<0):
        shift[rg:,:y-up] = star_rmbkg[:x-rg,up:]
    else:
        pass
    shift_part = np.zeros(shift.shape)
    shift_part[850:925,1050:1180] = shift[850:925,1050:1180]
    '''
    main_rmbkg = (data_main/exp_main - med)
    star_rmbkg[np.where(star_rmbkg<0.005)] = 0
    plt.imshow(main_rmbkg,vmin=0.005,vmax=0.5,origin='center')
    plt.xlim(1100,1250)
    plt.ylim(850,1000)
    plt.show()
    diff = main_rmbkg-shift
    plt.imshow(diff,vmin=0,vmax=0.02,origin='center')
    plt.xlim(1100,1250)
    plt.ylim(850,1000)
    plt.show()
    '''
    main_rmstar = (data_main/exp_main - shift)
    if fill == True:
        main_rmstar = (main_rmstar - 1000*shift_part)
    #main_rmstar[np.where(main_rmstar<0)] = med[np.where(main_rmstar<0)]
    # shift the main image
        angle = 35.1536*np.pi/180
        patch_d = 15
        patch = np.zeros(data_main.shape)
        x, y = data_main.shape
        p_rg = round(patch_d*np.cos(angle))
        p_up = round(patch_d*np.sin(angle))
        patch[p_rg:,p_up:] = data_main[:x-p_rg,:y-p_up]
        patch = patch/exp_main
        y_star,x_star = np.where(shift_part!=0)
        main_rmstar[np.where(shift_part!=0)] = patch[np.where(shift_part!=0)]
    ###
        if fill_sed == True:
            patch2 = np.zeros(main_rmstar.shape)
            x, y = main_rmstar.shape
            patch2[p_rg:,p_up:] = main_rmstar[:x-p_rg,:y-p_up]
            main_rmstar[np.where(shift_part!=0)] = patch2[np.where(shift_part!=0)]
    ###
    #plt.imshow(data_main, vmin=0, vmax=0.1*exp_star, origin='bottom')
    #plt.show()
    #plt.imshow(main_rmstar, vmin=0, vmax=0.1, origin='bottom')
    #plt.plot(x_star, y_star, 'rx',MarkerSize=1)
    #plt.show()
    #if fill == True:
    #    plt.imshow(patch, vmin=0, vmax=0.1, origin='bottom')
    #    plt.plot(x_star, y_star, 'rx',MarkerSize=1)
    #    plt.show()
    main_rmstar = main_rmstar*exp_main
    return main_rmstar

def read_motion():
    file_name = 'motion_log.txt'
    file_path = get_path('../docs/'+file_name)
    data = pd.read_csv(file_path, sep=' ', header=0)
    #data_group = data.groupby('ASTEROID')
    motion_dict = {}
    for i in range(0,len(data)):
        obsid = '000'+str(data['OBS_ID'][i])
        aster = data['ASTEROID'][i]
        motion_dict[obsid] = {'aster':aster,
                              'angle':data['ANGLE'][i],
                              'motion_p':data['MOTION_P'][i],
                              'motion_v':data['MOTION_V'][i]}
    return motion_dict

motion_dict = read_motion()

def rotateImage(img, angle, pivot):
    import scipy.ndimage as ndimage
    padX = [img.shape[1] - pivot[1], pivot[1]]
    padY = [img.shape[0] - pivot[0], pivot[0]]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, angle, reshape=False, order = 1,mode = 'constant',cval = 0.)
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]

def flipImage(img, pivot):
    import scipy.ndimage as ndimage
    padX = [img.shape[1] - pivot[1], pivot[1]]
    padY = [img.shape[0] - pivot[0], pivot[0]]
    imgP = np.pad(img, [padY, padX], 'constant')
    #imgR = ndimage.rotate(imgP, angle, reshape=False, order = 1,mode = 'constant',cval = 0.)
    imgF = np.flipud(imgP)
    return imgF[padY[0] : -padY[1], padX[0] : -padX[1]]

def replace(aster, obsid, regs, shift_v=0, motion_dict=motion_dict,radec='fitting'):
    # obtain flipped image
    print(aster,obsid)
    hdu = read_data(aster, obsid)
    data = hdu[1].data
    exp = hdu[1].header['EXPOSURE']
    if radec=='fitting':
        ra,dec=anker2sky(aster, obsid, 'colina96', '_smeargauss', output='radec')
    elif radec=='header':
        ra = hdu[1].header['RA_OBJ']
        dec = hdu[1].header['DEC_OBJ']
    else:
        pass
    wS =wcs.WCS(header=hdu[1].header,key='S',relax=True,)
    ypix, xpix = wS.wcs_world2pix(ra,dec,0)
    xpix = round(float(xpix))
    ypix = round(float(ypix))
    angle = motion_dict[obsid]['angle']
    data_rota1 = rotateImage(data,-angle,(xpix,ypix))
    # flip and shift
    xpix_flip = xpix+shift_v
    #plt.imshow(data_rota1, vmin=0, vmax=0.1*exp, origin='bottom')
    #plt.plot(800,xpix_flip,'kx',MarkerSize=5)
    #plt.xlim(600,1000)
    #plt.ylim(600,700)
    #plt.show()
    data_flip = flipImage(data_rota1,(xpix_flip,ypix))
    #plt.imshow(data_flip, vmin=0, vmax=0.1*exp, origin='bottom')
    #plt.plot(800,xpix_flip,'kx',MarkerSize=5)
    #plt.xlim(600,1000)
    #plt.ylim(600,700)
    #plt.show()
    data_fill = np.copy(data)
    med  = bkg_med(aster)
    motion_p = motion_dict[obsid]['motion_p']
    motion_v = motion_dict[obsid]['motion_v']
    shift_p = abs(round(motion_p))
    #print(shift_p)
    #plt.imshow(data, vmin=0, vmax=0.1*exp, origin='bottom')
    #plt.show()
    for reg in regs:
        if reg[4] == 'up':
            posi = 1
        elif reg[4] == 'dw':
            posi = -1
        elif (reg[4] == 'cn') or (reg[4] == 'bk'):
            posi = 0
        dirc = (motion_p/abs(motion_p))*(motion_v/abs(motion_v))*posi
        data_flip_shiftp = np.zeros(data_flip.shape)
        if dirc == 1:
            data_flip_shiftp[:,shift_p:] = data_flip[:,:data_flip.shape[1]-shift_p]
        elif dirc == -1:
            data_flip_shiftp[:,:data_flip.shape[1]-shift_p] = data_flip[:,shift_p:]
        elif dirc == 0:
            data_flip_shiftp = data_flip
        # rotate
        data_rota2 = rotateImage(data_flip_shiftp,angle,(xpix,ypix))
        plt.imshow(data_rota2, vmin=0, vmax=0.1*exp, origin='bottom')
        plt.show()
        # identify stars
        star = data - np.median([med,data,data_rota2],axis=0)
        star_part = np.zeros(star.shape)
        star_part[reg[0]:reg[1],reg[2]:reg[3]] = star[reg[0]:reg[1],reg[2]:reg[3]]
        star_part[np.where(star_part<0.005*exp)] = 0
        #plt.imshow(star_part, vmin=0, vmax=0.1*exp, origin='bottom')
        #plt.show()
        # fill
        if reg[4] == 'bk':
            data_fill[np.where(star_part!=0)] = med[np.where(star_part!=0)]*exp
        else:
            data_fill[np.where(star_part!=0)] = data_rota2[np.where(star_part!=0)]
    #
    #plt.imshow(data_fill, vmin=0, vmax=0.1*exp, origin='bottom')
    #plt.show()
    return data_fill

def replace_img(obsid):
    repl_dict = {'00091525001':{'aster':'iris','regs':[[900,980,1000,1080,'up']],'radec':'fitting','shift_v':0},
                 '00091519001':{'aster':'iris','regs':[[960,1000,1020,1080,'up'],[810,850,1230,1270,'up']],'radec':'fitting','shift_v':0},
                 '00091523001':{'aster':'iris','regs':[[930,990,1090,1150,'up']],'radec':'fitting','shift_v':0},
                 '00091532001':{'aster':'nysa','regs':[[1050,1090,860,900,'up'],[835,860,1130,1180,'dw']],'radec':'fitting','shift_v':-5},
                 '00091593001':{'aster':'themis','regs':[[900,930,1080,1120,'up']],'radec':'header','shift_v':0},
                 '00091595001':{'aster':'themis','regs':[[1050,1090,900,940,'up']],'radec':'header','shift_v':0},
                 '00091027001':{'aster':'dembowska','regs':[[970,1000,920,960,'cn'],[990,1010,880,910,'dw'],[820,860,1130,1180,'cn'],[760,800,1240,1290,'up']],'radec':'fitting','shift_v':0},
                 '00091505001':{'aster':'flora','regs':[[720,760,1150,1200,'dw']],'radec':'fitting','shift_v':-1},
                 #'00091503001':{'aster':'flora','regs':[[820,850,1050,1090,'cn']],'radec':'header','shift_v':3},
                 '00091501001':{'aster':'flora','regs':[[910,980,880,960,'cn'],[830,870,1030,1080,'cn']],'radec':'fitting','shift_v':-2},
                 '00091559001':{'aster':'hygiea','regs':[[880,910,1090,1120,'dw'],[880,920,1130,1190,'up']],'radec':'header','shift_v':2},
                 '00091521001':{'aster':'iris','regs':
                                [[1010,1050,1000,1040,'up'],[970,1020,1030,1090,'up'],
                                 [930,980,1070,1180,'up'],[900,940,1120,1220,'up'],
                                 [860,910,1170,1280,'up'],[820,870,1224,1324,'up'],
                                 [770,830,1280,1390,'up'],[724,780,1350,1454,'up'],
                                 [950,990,1010,1040,'dw']],'radec':'fitting','shift_v':-3},
                }
    subt_dict = {'00091268001':{'aster':'flora','star':'00091501001','offset':(5,5)},
                 '00091507001':{'aster':'flora','star':'00091505001','offset':(2,4)},
                 '00091503001':{'aster':'flora','star':'00091501001','offset':(-5,-7)},
                 '00091207002':{'aster':'juno','star':'00091206002','offset':(-2,1)},
                }
    if obsid in repl_dict.keys():
        inp_dict = repl_dict[obsid]
        data_repl = replace(inp_dict['aster'],obsid,inp_dict['regs'],radec=inp_dict['radec'],shift_v=inp_dict['shift_v'])
        return data_repl
    elif obsid in subt_dict.keys():
        inp_dict = subt_dict[obsid]
        data_repl = shift_star(inp_dict['aster'],obsid,inp_dict['star'],offset=inp_dict['offset'])
        return data_repl
    else:
        return None

#replace('iris','00091525001',regs=[[900,980,1000,1080,'up']]) #----
#replace('iris', '00091519001',regs=[[960,1000,1020,1080,'up'],[810,850,1230,1270,'up']])#,shift_v=-2) ----
#replace('iris', '00091523001',regs=[[930,990,1090,1150,'up']]) #----
#replace('nysa', '00091532001',regs=[[1050,1090,860,900,'up'],[835,860,1130,1180,'dw']],shift_v=-5) #----

replace('themis', '00091593001',regs=[[900,930,1080,1120,'up']],radec='header') #----
#replace('themis', '00091595001',regs=[[1050,1090,900,940,'up']],radec='header') #----
#replace('vesta', '00091022002',regs=[[850,880,1140,1170,'up']],radec='header',shift_v=1) #----
#replace('dembowska', '00091027001',regs=[[970,1000,920,960,'cn'],[990,1010,880,910,'dw'],[820,860,1130,1180,'cn'],[760,800,1240,1290,'up']]) #----
#x replace('dembowska', '00091237001',regs=[[850,910,1070,1170,'cn']],shift_v=-4) # something left
#replace('flora', '00091505001',regs=[[720,760,1150,1200,'dw']],shift_v=-1) #----
#x replace('flora', '00091268001',regs=[[790,860,1040,1100,'up']],shift_v=-3) #x center
#x replace('flora', '00091507001',regs=[[710,760,1200,1250,'cn']],radec='header',shift_v=-1) # something left
#x replace('flora', '00091503001',regs=[[820,850,1050,1090,'cn']],radec='header',shift_v=3) # x center
#replace('iris', '00091521001',regs=[[1010,1050,1000,1040,'up'],[970,1020,1030,1090,'up'],
#                                    [930,980,1070,1180,'up'],[900,940,1120,1220,'up'],
#                                    [860,910,1170,1280,'up'],[820,870,1224,1324,'up'],
#                                    [770,830,1280,1390,'up'],[724,780,1350,1454,'up'],
#                                    [950,990,1010,1040,'dw']],shift_v=-3) #----
#x replace('juno', '00091207002',regs=[[780,850,1150,1220,'up']],radec='header',shift_v=1) # x center
#replace('hygiea', '00091559001',regs=[[880,910,1090,1120,'dw'],[880,920,1130,1190,'up']],radec='header',shift_v=2) #----
#x replace('nysa', '00091534001',regs=[[900,980,1000,1080,'up']]) #x center
#x replace('nysa', '00091540001',regs=[[900,980,1000,1080,'up']]) #x center
#x replace('nysa', '00091538001',regs=[[900,980,1000,1080,'up']]) #x center
#replace('flora', '00091501001',regs=[[910,980,880,960,'cn'],[830,870,1030,1080,'cn']],shift_v=-2)

#shift_star('hygiea', '00091556001', '00091554001')
#shift_star('massalia', '00091545001', '00091543001')
#shift_star('dembowska', '00091237001', '00091239001')

#shift_star('ceres','00091242001','00091244001') #s
#shift_star('dembowska','00091027001','00091237001') #s '00091237001','00091239001'
#shift_star('dembowska','00091237001','00091027001') #l '00091027001','00091239001','00091241001'


# hygiea: 00091552001, 00091556001, 00091559001, 00091554001, 00091550001
# iris: 00091527001, 00091523001, 00091519001, 00091521001, 00091525001
# nysa: 00091532001, 00091536001, 00091538001, 00091540001, 00091534001
# themis: 00091595001, 00091591001, 00091581001, 00091589001, 00091593001, 00091583001
# vesta: 00091022002, 00091208002, 00091198002, 00091197002
# dembowska: 00091239001, 00091237001, 00091241001, 00091027001
# flora: 00091503001, 00091507001, 00091268001, 00091505001, 00091501001
# juno: 00091026003, 00091204002, 00091205002, 00091206002, 00091203002, 00091207002
# edge
#shift_star('hygiea', '00091559001', '00091556001') # no other obs has similar pointing, or rotation?
#shift_star('iris', '00091525001', '00091525001') #l rotation, edge/out
#shift_star('iris', '00091519001', '00091519001') #rotation, edge/out
#shift_star('iris', '00091523001', '00091523001') #rotation, out
#shift_star('iris', '00091527001', '00091527001') #rotation, out, !
#shift_star('nysa', '00091532001', '00091532001') #rotation, out
#shift_star('themis', '00091593001', '00091593001') #rotation, out
#shift_star('themis', '00091595001', '00091595001') #rotation, out
#shift_star('vesta', '00091022002', '00091198002') #little offset x

# whole
#shift_star('dembowska', '00091027001', '00091237001') #s little offset
#shift_star('dembowska', '00091237001', '00091239001') #l little offset
#shift_star('flora', '00091505001', '00091507001') #s little offset x
#shift_star('flora', '00091268001', '00091501001',offset=(5,5)) #l little offset
#shift_star('flora', '00091507001', '00091505001',offset=(2,4)) #l little offset
#shift_star('flora', '00091503001', '00091501001', offset=(-5,-7)) #m little offset
#shift_star('iris', '00091521001', '00091521001') #rotation, edge/out,line
#shift_star('juno', '00091207002', '00091206002', offset=(-2,1)) # little offset
#shift_star('nysa', '00091534001', '00091534001') #rotation, center
#shift_star('nysa', '00091540001', '00091540001') #center, edge, line
#shift_star('nysa', '00091538001', '00091538001') #center
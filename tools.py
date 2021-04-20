from _mypath import thisdir
from os import listdir
from scipy import interpolate
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astropy import units as u
from astropy.coordinates import SkyCoord
import os

def get_path(*rela_path, to_file=False):
    """get an absolute path
    if to_file=True, the path will link to a .img.gz file according to inputted OBS_ID and FILTER
    if to_file=False, the path will link to any file/dir according to the rela_path
    
    Inputs: 
    1> 'OBS_ID', 'FILTER'('uw1', 'uw2', 'uvv') 
    2> 'rela_path', to_file=F
    Outpus: 
    absolute path, string
    """
    if to_file == True:
        data_path = '../data/'
        data_abs_path = os.path.join(thisdir, data_path)
        folders = [one for one in listdir(data_abs_path) if 'raw' in one]
        for folder in folders:
            folder_abs_path = ''.join([data_abs_path, folder, '/'])
            obs_list = listdir(folder_abs_path)
            obs_list = [obs for obs in obs_list if obs[0] != '.']
            if rela_path[0] in obs_list:
                break
        map_path = '/uvot/image/'
        file_name = 'sw' + rela_path[0] + rela_path[1] + '_sk.img.gz'
        file_path = ''.join([folder_abs_path, rela_path[0], map_path, file_name])
        return file_path
    else:
        return os.path.join(thisdir, rela_path[0])

def error_prop(method, x, x_err, y, y_err):
    if method == 'sum' or method == 'sub':
        return np.sqrt(np.power(x_err,2)+np.power(y_err,2))
    elif method == 'mul':
        err_1 = np.power(x_err*y,2)
        err_2 = np.power(y_err*x,2)
        return np.sqrt(err_1+err_2)
    elif method == 'div':
        err_1 = np.power(x_err/y,2)
        err_2 = np.power(x*y_err/np.power(y,2),2)
        return np.sqrt(err_1+err_2)

def as2au(arcsec, dis_au):
    return dis_au*2*np.pi*arcsec/(3600.*360)
    
def au2km(au):
    return 149597870.7*au

def integ(x, y, step_num=False):
    if step_num:
        f = interpolate.interp1d(x, y, fill_value='extrapolate')
        x = np.linspace(x[0], x[-1], num=step_num)
        y = f(x)
    result = 0
    for i in range(len(x)):
        result += x[i]*y[i]
    return result

def getRadec(obsid,horizon_id,coord,aster):
    # coord: '@swift' or 500
    img_name = 'sw'+str(obsid)+'ugu_sk.img.gz'
    #img_name = 'sw'+obsid+'uw2_sk.img'
    img_path = '/Users/zexixing/Research/swiftASTER/data/'+aster+'/'+\
               obsid+'/uvot/image/'+img_name
    hdul = fits.open(img_path)
    start = Time(hdul[0].header['DATE-OBS'])
    end = Time(hdul[0].header['DATE-END'])
    dt = end-start
    mid_time = start + 1/2*dt
    obj = Horizons(id=horizon_id,
                   location=coord,
                   epochs=mid_time.jd)
    eph = obj.ephemerides()[0]
    ra = eph['RA']
    dec = eph['DEC']
    return ra, dec

def hmsdms2deg(hmsdms):
    c = SkyCoord(hmsdms, unit=(u.hourangle, u.deg)) # '00 42 30 +41 12 00'
    return c.ra.value, c.dec.value

print(au2km(1))
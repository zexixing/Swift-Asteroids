import numpy as np
from astropy.io import fits
from itertools import islice, product
from tools import get_path
from conversion import *
import math

def load_img(img_name,relative_path='', ext=0):
    '''
    load fits data
    inputs：name of a fits file in '../docs/', str
    return: array
    '''
    img_path = get_path('../docs/'+relative_path+img_name)
    if isinstance(ext,int)==True:
        img_data = fits.open(img_path)[ext].data
        return img_data
    else:
        dic_data = {}
        fits_file = fits.open(img_path)
        # ext = {0:'img',1:'err2',2:'exp'}
        for i in ext:
            dic_data[ext[i]] = fits_file[i].data
        return dic_data

def load_header(img_name,relative_path=''):
    '''
    load fits header
    inputs: str
    return: dictionary
    '''
    img_path = get_path('../docs/'+relative_path+img_name)
    img_header = fits.open(img_path)[1].header #~~TODO:
    return img_header

def get_dist(point_1, point_2):
    '''
    get distance between two points
    inputs: pixel position from ds9
    reutrn: float
    '''
    i2 = (point_1[0]-point_2[0])**2
    j2 = (point_1[1]-point_2[1])**2
    dist = np.sqrt(i2+j2)
    return dist

def load_reg_list(reg_name,relative_path=''):
    '''
    load regions from a .reg file (requirement: a list of circle regions)
    inputs: .reg file name (str)
    outputs: y position list; x position list; radius list
    '''
    reg_path = get_path('../docs/'+relative_path+reg_name)
    reg_file = open(reg_path)
    center_i_list = []
    center_j_list = []
    radiu_list = []
    for line in islice(reg_file, 3, None):
        reg_data = line.split('(')[1]
        reg_data = reg_data.split(')')[0]
        reg_data = reg_data.split(',')
        center_j_list.append(float(reg_data[0]))
        center_i_list.append(float(reg_data[1]))
        if len(reg_data[2:]) == 1:
            radiu_list.append(float(reg_data[2:][0]))
        else:
            radiu_list.append([float(k) for k in reg_data[2:]])
    return center_i_list, center_j_list, radiu_list

def mask(img_name, lowlim, uplim, relative_path=''):
    '''
    #not used in final pipeline
    mask pixels with values larger than uplim and lower than lowlim
    inputs: name of the image to be masked (str); float; float
    return: mask array
    '''
    # 1: unmasked, 0: masked
    img_data = load_img(img_name,relative_path)
    mask_img = np.ones(img_data.shape)
    i_range = np.arange(len(img_data))
    j_range = np.arange(len(img_data[0]))
    for i in i_range:
        for j in j_range:
            if img_data[i,j]>uplim or img_data[i,j]<lowlim:
                mask_img[i,j] = 0
    return mask_img

def limit_loop(img_name, center, r, relative_path='',ext=0):
    '''
    crop the image to get a square with sides of 2r
    to reduce loop times in the following aperture photometry codes
    inputs: image name (str) or image array data (array); 
            the center of the region to be cut out (values from ds9)((x,y))
            length of side: 2r (float)
    return: the whole data array (array); 
            the center of the region to be cut out (to be used in python)((x,y))
            x position list of the cropped data array (list)
            y position list of the cropped data array (list)
    '''
    if isinstance(img_name, str) == True:
        img_data = load_img(img_name,relative_path,ext=ext)
    else:
        img_data = img_name
    if isinstance(ext,int)==True:
        i_range = np.arange(len(img_data))
        j_range = np.arange(len(img_data[0]))
    else:
        i_range = np.arange(len(img_data['img']))
        j_range = np.arange(len(img_data['img'][0])) 
    cen_pix = [center[1]-1, center[0]-1]
    i_range = i_range[i_range >= cen_pix[0]-r]
    i_range = i_range[i_range <= cen_pix[0]+r]
    j_range = j_range[j_range >= cen_pix[1]-r]
    j_range = j_range[j_range <= cen_pix[1]+r]
    return img_data, cen_pix, i_range, j_range

def mask_region(mask_img, mask_name, relative_path=''):
    if isinstance(mask_img, str):
        img_data = load_img(mask_img,relative_path)
        mask_img = np.ones(img_data.shape)
    regions = load_reg_list(mask_name,relative_path)
    x_list = regions[1]
    y_list = regions[0]
    r_list = regions[2]
    center_list = list(zip(x_list, y_list))
    for k in range(len(r_list)):
        img_data, cen_pix, i_range, j_range = limit_loop(mask_img, center_list[k], r_list[k], relative_path)
        for i in i_range:
            for j in j_range:
                pos_pix = [i, j]
                pos_cen_dist = get_dist(cen_pix, pos_pix)
                if pos_cen_dist < r_list[k]:
                    mask_img[i,j] = 0
    return mask_img


# new modify
def phot_ct(img_name, center, src_r, start=0, step = 2, shape='total',
            method='mean', mask_img=False, relative_path='',ext=0):
    # aper list
    step_num = math.ceil((src_r-start)/step)
    aper_list = np.linspace(start,step*step_num,step_num+1)
    aper_list = (aper_list[1:]+aper_list[:-1])/2.
    if start == False:
        start = 0
    img_data, cen_pix, i_range, j_range = \
        limit_loop(img_name, center, src_r, relative_path, ext)
    if isinstance(ext,int) == True:
        exp_data = np.ones(img_data.shape)
    else:
        dic_data = img_data
        data_keys = list(dic_data.keys())
        img_data = dic_data['img']
        if 'exp' in data_keys:
            exp_data = dic_data['exp']
        else:
            exp_data = np.ones(img_data.shape)
    if isinstance(mask_img, bool):
        if mask_img == False:
            mask_img = np.ones(img_data.shape)
    pos_data = np.argwhere(img_data!=np.nan)
    dist_data = np.sqrt(np.sum(np.power(pos_data-cen_pix,2),axis=1)).reshape(img_data.shape) #TODO:
    index_data = np.ceil((np.array(dist_data)-start)/step)-1
    index_data[np.where(index_data==-1)]=0
    # mask
    import numpy.ma as ma
    mask_map = ma.masked_where(dist_data>src_r,img_data)
    if start !=0:
        mask_map = ma.masked_where((dist_data<=start),mask_map)
    else:
        pass
    mask_map_unmask = mask_map
    mask_map = ma.masked_where(exp_data!=np.max(exp_data),mask_map)
    if isinstance(mask_img, str):
        if mask_img == 'star':
            mask_map = ma.masked_where(img_data>1.5,mask_map)
            mask_img = np.ones(img_data.shape)
    mask_map = ma.masked_where(mask_img==0,mask_map)
    #plt.imshow(mask_map,vmin=0,vmax=3)
    #plt.show()
    # get median profile
    image_list = img_data[mask_map.mask==False]
    index_list = index_data[mask_map.mask==False]
    index_unmask_list = index_data[mask_map_unmask.mask==False]
    # judge aperture
    if len(list(set(index_data.flatten())))<len(aper_list):
        aper_list = aper_list[:len(list(set(index_data.flatten())))]
        print('Given aperture is so large that a masked/blank annuli is included.')
        print('The area has been cut.')
    # get results
    median_list = []
    pixel_mask_list = []
    pixel_unmask_list = []
    err_list = []
    if shape == 'total':
        pixel_unmask = len(mask_map_unmask.mask[mask_map_unmask.mask==False])
        pixel_mask = len(mask_map.mask[mask_map.mask==False])
        if method == 'mean':
            count = np.sum(image_list)/pixel_mask*pixel_unmask
            err = np.sqrt(np.sum(image_list))/pixel_mask*pixel_unmask
            return count, pixel_unmask, pixel_mask, err
        elif method == 'median':
            count = np.median(image_list)*pixel_unmask
            err = np.std(image_list)*pixel_unmask
            return count, pixel_unmask, pixel_mask, err
    elif shape == 'azim':
        for ind in np.arange(np.min(index_list),np.max(index_list)+1):
            median_unit = image_list[index_list==ind]
            median_list.append(np.median(median_unit))
            err_list.append(np.std(median_unit))
            pixel_mask_list.append(len(median_unit))
            pixel_unmask_list.append(len(index_unmask_list[index_unmask_list==ind]))
        median_list = np.array(median_list)
        err_list = np.array(err_list)
        pixel_unmask_list = np.array(pixel_unmask_list)
        pixel_mask_list = np.array(pixel_mask_list)
        if len(pixel_unmask_list) == len(pixel_mask_list):
            if method == 'median':
                count = np.sum(median_list*pixel_unmask_list)
                err = np.sum(err_list*pixel_unmask_list)
                return count, np.sum(pixel_unmask_list), np.sum(pixel_mask_list), err
            elif method == 'mean':
                count = np.sum(sum_list/pixel_mask_list*pixel_unmask_list)
                err = np.sum(np.sqrt(sum_list)/pixel_mask_list*pixel_unmask_list)
                return count, np.sum(pixel_unmask_list), np.sum(pixel_mask_list), err
        else:
            raise Exception('some index was removed when masking img!')
    else:
        raise Exception('Invalid input of shape!')

def multi_circle_ct(img_name, center_list, r_list, method='mean', mask_img = False, relative_path=''):
    count_list = []
    pixel_list = []
    for n in range(len(r_list)):
        result = phot_ct(img_name, center_list[n], r_list[n], method=method, mask_img=mask_img, relative_path=relative_path)
        count_list.append(result[0])
        pixel_list.append(result[1])
    return count_list, pixel_list

def load_reg_list(reg_name,relative_path=''):
    reg_path = get_path('../docs/'+relative_path+reg_name)
    reg_file = open(reg_path)
    center_i_list = []
    center_j_list = []
    radiu_list = []
    for line in islice(reg_file, 3, None):
        reg_data = line.split('(')[1]
        reg_data = reg_data.split(')')[0]
        reg_data = reg_data.split(',')
        center_j_list.append(float(reg_data[0]))
        center_i_list.append(float(reg_data[1]))
        if len(reg_data[2:]) == 1:
            radiu_list.append(float(reg_data[2:][0]))
        else:
            radiu_list.append([float(k) for k in reg_data[2:]])
    return center_i_list, center_j_list, radiu_list

def reg2bg_bri(img_name, bg_method, bg_center, bg_r, count_method, mask_img = False, relative_path='', ext=0):
    if bg_method == 'single':
        result = phot_ct(img_name, bg_center, bg_r, method=count_method, mask_img=mask_img, relative_path=relative_path, ext=ext)
    elif bg_method == 'donut':
        result = phot_ct(img_name, bg_center, bg_r[1], start=bg_r[0], method=count_method, mask_img=mask_img, relative_path=relative_path, ext=ext)
    elif bg_method == 'multi':
        result = multi_circle_ct(img_name, bg_center, bg_r, count_method, mask_img, relative_path=relative_path, ext=ext)[:2]
        bg_count = np.array(result[0])
        bg_pixel = np.array(result[1])
    else:
        raise Exception('check the input of method (single/dount/multi)')
    if count_method == 'mean':
        if bg_method == 'multi':
            bg_bri = bg_count/bg_pixel
            bg_bri = np.mean(bg_bri)
            n = len(bg_center)
            bg_bri_err = (1/n)*np.sqrt(np.sum(abs(bg_count)/np.power(bg_pixel, 2)))
        else:
            bg_bri = result[0]/result[1]
            bg_bri_err = np.sqrt(abs(result[0]))/result[1]
    elif count_method == 'median':
        if bg_method == 'multi':
            bg_bri = bg_count/bg_pixel
            bg_bri = np.median(bg_bri)
            bg_bri_err = np.std(bg_bri)
        else:
            bg_bri = result[0]/result[1]
            bg_bri_err = result[2]/result[1]#~~TODO:
    else:
        raise Exception('check the input of method (mean/median)')
    return bg_bri, bg_bri_err

def aper_phot_cr(src_count, src_pixel, 
                 bg_bri, exposure, 
                 src_count_err, bg_bri_err):
    """aperture photometry
    """
    bg_count = bg_bri*src_pixel
    bg_count_err = bg_bri_err*src_pixel
    net_count = src_count - bg_count
    net_count_err = np.sqrt(bg_count_err**2
                            +src_count_err**2)
    net_cr = net_count/exposure
    net_cr_err = net_count_err/exposure
    snr = abs(net_cr/net_cr_err)
    return net_cr, net_cr_err, snr

def aper_phot(img_name, filt,
              src_center, src_r,
              bg_center, bg_r,
              src_method, bg_method,
              step=5, mask_img = False, 
              start = False, relative_path='',ext=0):
    '''
    src_method = 'total_mean' or 'total_median'
                 or 'azim_mean' or 'azim_median'
    bg_method = 'single_mean' or 'single_median'
                or 'donut_mean' or 'donut_median'
                or 'multi_mean'
    
    return cr, cr_err, snr, mag, mag_err, bg_cr, bg_cr_err
    '''
    # src photometry: get src_count, src_pixel
    src_shape = src_method.split('_')[0]
    src_stat = src_method.split('_')[1]
    result = phot_ct(img_name,
                     src_center, src_r, start,
                     step = step,
                     shape = src_shape, 
                     method = src_stat, 
                     mask_img = mask_img, 
                     relative_path = relative_path,
                     ext = ext)
    src_count = result[0]
    src_pixel_unmask = result[1]
    src_pixel_mask = result[2]
    src_err = result[3]
    #print('src_count:'+str(src_count)+'; src_pixel_unmask:'+str(src_pixel_unmask))
    # bg photometry: get bg_bri
    bg_shape = bg_method.split('_')[0]
    bg_stat = bg_method.split('_')[1]
    bg_bri, bg_bri_err = reg2bg_bri(img_name, bg_shape, 
                                    bg_center, bg_r, 
                                    bg_stat, mask_img = mask_img, 
                                    relative_path=relative_path,
                                    ext=ext)
    #print('bg_bri:'+str(bg_bri))
    # photometry
    if filt:
        exposure = float(load_header(img_name,relative_path)['EXPOSURE'])#EXPTIME
    else:
        exposure = 1.
    cr, cr_err, snr = aper_phot_cr(src_count, src_pixel_unmask,
                                   bg_bri, exposure,
                                   src_err, bg_bri_err)
    if filt:
        mag, mag_err = cr2mag(cr, cr_err, filt)
    else:
        mag = float('NaN')
        mag_err = float('NaN')
    bg_cr = bg_bri/exposure
    bg_cr_err = bg_bri_err/exposure
    return (cr, cr_err), snr, (mag, mag_err), (bg_cr, bg_cr_err)


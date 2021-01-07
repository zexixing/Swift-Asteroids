import numpy as np

def getCalData(wheelpos,chatter=3,s=0,calfile=None,caldir=None, msg=''):
   '''Retrieve the calibration data for the anchor and dispersion (wavelengths).
   
   Parameters
   ----------
   Xphi, Yphi : float
      input angles in degrees, from, e.g., `findInputAngle`.
      
   wheelpos : int, {160,200,955,1000}
      filter wheel position selects grism
       
   date : swifttime in seconds
      obsolete - not used
      
   kwargs : dict
     optional arguments
      
     - **calfile** : str
     
       calibration file name
     
     - **caldir** : str
     
       path of directory calibration files
     
     - **mode** : str
     
       interpolation method. Use 'bilinear' only.
     
     - **kx**, **ky** : int, {1,2,3}
       order of interpolation. Use linear interpolation only. 
       
     - **s** : float
       smoothing factor, use s=0.
       
     - **chatter** : int
       verbosity
   
   Returns
   -------
   anker, anker2 : list
      coordinate of anchor in first order.
       
   C_1, C_2 : 
      dispersion in first and second order.
       
   theta : float
      find angle of dispersion on detector as 180-theta.
      
   data : FITS_rec
      the wavecal data table  
   
   Notes
   -----
   Given the input angle Xphi, Yphi in deg., the filterwheel 
   position, and the date the spectrum was taken (in swift seconds), 
   this gets the calibration data. 
      
   The boresight must be set to the one used in deriving the calibration. 
      
   '''

   import os
   import numpy as np
   try:
      from astropy.io import fits as pyfits
   except:   
      import pyfits
   from scipy import interpolate
   
   #==================================================================
   # The following calculation in reverse prepared the zemax model for 
   # the calibration table lookup. Keep for the record. UV Nominal case.
   #  first calculate the offset of the rotate the input angles due to 
   #  the difference in boresight of grism and model  
   #   = input_angle + grism_bs_angle - model_bs
   # scale = 6554.0  (deg/pix)
   # xfi = Xphi + (( 928.53-27) - (1100.5+8))/scale
   # yfi = Yphi + ((1002.69- 1) - (1100.5-4))/scale
   # rx,ry = uvotmisc.uvotrotvec(xf,yf,-64.6)
   #==================================================================
   if calfile == None:
      #
      #  get the calibration file
      #
      try:
         uvotpy = os.getenv('UVOTPY')+'/uvotpy'
         caldb  = os.getenv('CALDB')
         if uvotpy != None: 
            caldir = uvotpy+'/calfiles'
         elif caldb != None:
            caldir = caldb+'/data/swift/uvota/bcf/grism/'
      except:
         print("CALDB nor UVOTPY environment variable set.")     
      
      #if caldir == None: 
      #   # hardcoded development system 
      #   caldir = '/Volumes/users/Users/kuin/dev/uvotpy.latest/calfiles'
         
      if wheelpos == 200: 
         calfile = 'swugu0200wcal20041120v001.fits'
         oldcalfile='swwavcal20090406_v1_mssl_ug200.fits'
         calfile = caldir+'/'+calfile
         if chatter > 1: print('reading UV Nominal calfile '+calfile)
      elif wheelpos == 160: 
         calfile='swugu0160wcal20041120v002.fits'
         oldcalfile= 'swwavcal20090626_v2_mssl_uc160_wlshift6.1.fits'
         calfile = caldir+'/'+calfile
         if chatter > 1: print('reading UV clocked calfile '+calfile) 
      elif wheelpos == 955: 
         calfile='swugv0955wcal20041120v001.fits'
         oldcalfile= 'swwavcal20100421_v0_mssl_vc955_wlshift-8.0.fits'
         calfile = caldir+'/'+calfile
         if chatter > 1: print('reading V Clockedcalfile '+calfile) 
      elif wheelpos == 1000: 
         calfile='swugv1000wcal20041120v001.fits'
         oldcalfile= 'swwavcal20100121_v0_mssl_vg1000.fits'
         calfile = caldir+'/'+calfile
         if chatter > 1: print('reading V Nominal calfile  '+calfile) 
      else:
          if chatter > 1: 
             print("Could not find a valid wave calibration file for wheelpos = ",wheelpos)
             print("Aborting")
             print("******************************************************************")
             raise IOError("missing calibration file")
               

   msg += "wavecal file : %s\n"%(calfile.split('/')[-1])
   #  look up the data corresponding to the (Xphi,Yphi) point in the 
   #  calibration file (which already has rotated input arrays) 
   #    
   cal = pyfits.open(calfile)
   if chatter > 0: print("opening the wavelength calibration file: %s"%(calfile))
   if chatter > 1: print(cal.info())
   hdr0 = cal[0].header
   hdr1 = cal[1].header
   data = cal[1].data
   # the rotated field grid xf,yf (inconsistent naming - use to be xrf,yrf)
   xf = xrf = data.field('PHI_X')
   N1 = int(np.sqrt( len(xf) ))
   if N1**2 != len(xf): 
      raise RuntimeError("GetCalData: calfile array not square" )
   if chatter > 2: print("GetCalData: input array size on detector is %i in x, %i in y"%(N1,N1))   
   xf = xrf = data.field('PHI_X').reshape(N1,N1)
   yf = yrf = data.field('PHI_Y').reshape(N1,N1)
   #  first order anchor and angle array
   xp1 = data.field('DETX1ANK').reshape(N1,N1)
   yp1 = data.field('DETY1ANK').reshape(N1,N1)
   return xf,yf,xp1,yp1

def bilinear(x1,x2,x1a,x2a,f,chatter=0):
   '''
   Given function f(i,j) given as a 2d array of function values at
   points x1a[i],x2a[j], derive the function value y=f(x1,x2) 
   by bilinear interpolation. 
   
   requirement: x1a[i] is increasing with i 
                x2a[j] is increasing with j
   20080303 NPMK                
   '''
   import numpy as np
   
   # check that the arrays are numpy arrays
   x1a = np.asarray(x1a)
   x2a = np.asarray(x2a)
      
   #  find the index for sorting the arrays
   n1 = len(x1a)
   n2 = len(x2a)
   x1a_ind = x1a.argsort()
   x2a_ind = x2a.argsort()
   
   #  make a sorted copy
   x1as = x1a.copy()[x1a_ind]
   x2as = x2a.copy()[x2a_ind]
   
   # find indices i,j for the square containing (x1, x2)
   k1s = x1as.searchsorted(x1)-1
   k2s = x2as.searchsorted(x2)-1
   
   #  find the indices of the four points in the original array
   ki = x1a_ind[k1s]
   kip1 = x1a_ind[k1s+1]
   kj = x2a_ind[k2s]
   kjp1 = x2a_ind[k2s+1]
   if chatter > 2:
       print('FIND solution in (x,y) = (',x1,x2,')')
       print('array x1a[k-5 .. k+5] ',x1a[ki-5:ki+5])
       print('array x2a[k-5 .. k+5] ',x2a[kj-5:kj+5])
       print('length x1a=',n1,'   x2a=',n2)
       print('indices in sorted arrays = (',k1s,',',k2s,')')
       print('indices in array x1a: ',ki, kip1)
       print('indices in array x2a: ',kj, kjp1)
      
   #  exception at border:
   if ((k1s+1 >= n1) ^ (k2s+1 >= n2) ^ (k1s < 0) ^ (k2s < 0) ):
      print('bilinear. point outside grid x - use nearest neighbor ')
      if ki + 1 > len(x1a) : ki = len(x1a) - 1
      if ki < 0 : ki = 0
      if kj + 1 > len(x2a) : kj = len(x2a) - 1
      if kj < 0 : kj = 0
      return f[ki, kj]
  
   # Find interpolated solution
   y1 = f[kj  ,ki  ]
   y2 = f[kj  ,kip1]
   y3 = f[kjp1,kip1]
   y4 = f[kjp1,ki  ]
    
   t = (x1 - x1a[ki])/(x1a[kip1]-x1a[ki])
   u = (x2 - x2a[kj])/(x2a[kjp1]-x2a[kj])
   
   y = (1.-t)*(1.-u)*y1 + t*(1.-u)*y2 + t*u*y3 + (1.-t)*u*y4
   if chatter > 2: 
      print('bilinear.                   x         y          f[x,y]    ')
      print('bilinear.   first  point ',x1a[ki  ],x2a[kj],  f[ki,kj])
      print('bilinear.   second point ',x1a[kip1],x2a[kj],  f[kip1,kj])
      print('bilinear.   third  point ',x1a[kip1],x2a[kjp1],  f[kip1,kjp1])
      print('bilinear.   fourth point ',x1a[ki  ],x2a[kjp1],  f[ki,kjp1])
      print('bilinear. fractions t, u ', t, u)
      print('bilinear. interpolate at ', x1, x2, y)
   return y    

def fileinfo(filestub,ext,lfilt1=None, directory='./',chatter=0, wheelpos=None, twait=40.0):
   '''finds files for spectrum, matching attitude and lenticular images 
      uncompresses gzipped files if found
      
      Parameters
      ----------
      filestub : str
        the base of the file name, with the Swift project convention,
        consisting of "sw" + the OBSID, i.e., "sw00032301001" 
      ext : int
        the number of the extension of the grism file to match
        
      kwargs : dict
      
       - **lfilt1** : str, optional
       
        name of first lenticular filter used. Must be one of 'uvw2',
        'uvm2','uvw1','u','b','v','wh'  
        
       - **directory** : path, str
       
        path for directory. This is the directory with the grism file.
        
       - **chatter** : int
       
        verbosity
                                
       - **twait** : float
       
        The maximum time allowed between begin and end of matched 
        exposures of grism-lenticular filter, for each match.
        
       - **wheelpos** : imt
        
        If given, use to discriminate between UV and Visual grisms.     
        
      Returns
      -------
      specfile, attfile: str
        filename of spectrum, 
        the attitude file name.
      
      lfilt1, lfilt2 : str
       lenticular filter file name before the grism exposure or None, 
       and the file name of the lenticular filter following the grism 
       
      lfilt1_ext,lfilt2_ext : int
       extension number for each lenticular filter matching exposure 

   '''
   import os
   try:
      from astropy.io import fits 
   except:
      import pyfits as fits
   from numpy import array
   
   ext_names =array(['uw2','um2','uw1','uuu','ubb','uvv','uwh'])
   lfiltnames=array(['uvw2','uvm2','uvw1','u','b','v','wh'])
   
   vvgrism = True
   uvgrism = True
   if wheelpos != None:
      if wheelpos < 500: 
         vvgrism = False
         specfile =  directory+filestub+'ugu_dt.img'
      else: 
         specfile =  directory+filestub+'ugv_dt.img'
         uvgrism = False
   
   if (not directory.endswith('/')) : 
      directory += '/' 
   auxildir = directory+'../../auxil/'
   attfile = None
      
   # test if u or v grism file and set variable 
   specfile = ' *filename not yet initialised (directory wrong?)* '
   if uvgrism & os.access(directory+filestub+'ugu_dt.img',os.F_OK):
        specfile =  directory+filestub+'ugu_dt.img'
        if chatter > 1: print('reading ',specfile)
   elif uvgrism & os.access(directory+filestub+'ugu_dt.img.gz',os.F_OK):
        specfile =  directory+filestub+'ugu_dt.img'
        os.system( 'gunzip '+specfile+'.gz' )
        if chatter > 1: print('reading ',specfile)
   elif vvgrism & os.access(directory+filestub+'ugv_dt.img',os.F_OK):     
        specfile =  directory+filestub+'ugv_dt.img'
        if chatter > 1: print('reading ',specfile)
   elif vvgrism & os.access(directory+filestub+'ugv_dt.img.gz',os.F_OK):     
        specfile =  directory+filestub+'ugv_dt.img'
        os.system( 'gunzip '+specfile+'.gz' )
        if chatter > 1: print('reading ',specfile)
   else:
        print("on call fileinfo(sw+obsid="+filestub+",ext=",ext,",lfilt1=",lfilt1,", directory="+directory,",wheelpos=",wheelpos,")")
        raise IOError("FILEINFO: cannot find %s: DET file not found - pls check directory/file provided  is correct" % specfile )
                
   #    attitude file:
   if os.access(directory+filestub+'pat.fits',os.F_OK):
        attfile =  directory+filestub+'pat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(directory+filestub+'pat.fits.gz',os.F_OK):
        attfile =  directory+filestub+'pat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   elif os.access(directory+filestub+'uat.fits',os.F_OK):     
        attfile =  directory+filestub+'uat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(directory+filestub+'uat.fits.gz',os.F_OK):     
        attfile =  directory+filestub+'uat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   elif os.access(directory+filestub+'sat.fits',os.F_OK):     
        attfile =  directory+filestub+'sat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(directory+filestub+'sat.fits.gz',os.F_OK):     
        attfile =  directory+filestub+'sat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'pat.fits',os.F_OK):
        attfile =  auxildir+filestub+'pat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'pat.fits.gz',os.F_OK):
        attfile =  auxildir+filestub+'pat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'uat.fits',os.F_OK):     
        attfile =  auxildir+filestub+'uat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'uat.fits.gz',os.F_OK):     
        attfile =  auxildir+filestub+'uat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'sat.fits',os.F_OK):     
        attfile =  auxildir+filestub+'sat.fits'
        if chatter > 1: print('found att file ',attfile)
   elif os.access(auxildir+filestub+'sat.fits.gz',os.F_OK):     
        attfile =  auxildir+filestub+'sat.fits'
        os.system( 'gunzip '+attfile+'.gz' )
        if chatter > 1: print('found att file ',attfile)
   #    filter file(s)
   lfilt1,lfilt2 = None,None 
   lfilt1_ext = None; lfilt2_ext=None
   hdu = fits.open(specfile)
   if len(hdu)-1 < ext: 
      raise IOError("Input error: extension not found in Grism file.")
   hdr = hdu[int(ext)].header   
   hdu.close()
   #hdr = fits.getheader(specfile,int(ext))
   tstart = hdr['TSTART']
   tstop  = hdr['TSTOP'] 
   if chatter > 1: 
      print('grism times : %s - %s '%(tstart,tstop))
   lfile=None
   #  
   for k in range(len(ext_names)):
      ftyp = ext_names[k]
      lfiltyp = lfiltnames[k]
      if chatter > 1: print("testting for "+directory+filestub+ftyp+'_sk.img')
      if os.access(directory+filestub+ftyp+'_sk.img',os.F_OK):
        lfile =  directory+filestub+ftyp+'_sk.img'
        if chatter > 1: 
           print('found lenticular sky file ',lfile) 
      elif os.access(directory+filestub+ftyp+'_sk.img.gz',os.F_OK):
        lfile =  directory+filestub+ftyp+'_sk.img' 
        os.system( 'gunzip '+lfile+'.gz' )
        if chatter > 1: print('found lenticular sky file ',lfile)
      if lfile != None: 
         # check if it has an extension before or after specfile[ext] 
         xxx = fits.open(lfile)
         for i in range(1,len(xxx)):
            t1 = xxx[i].header['TSTART']
            t2 = xxx[i].header['TSTOP']
            if abs(t2-tstart) < twait:
               lfilt1 = lfiltyp
               lfilt1_ext = i
               if chatter > 0: print("lenticular file observation preceeding grism observation")
            if abs(t1-tstop) < twait:
               lfilt2 = lfiltyp
               lfilt2_ext = i   
               if chatter > 1: print("lenticular file observation after grism observation")
         lfile = None 
         xxx.close()
              
   # wrapup in case there is only one, but in lfilt2.
   if ((lfilt1 == None) & (lfilt2 != None)): 
      if chatter > 2: print("putting only filter info in filter 1")
      lfilt1 = lfilt2
      lfilt2 = None
      lfilt1_ext = lfilt2_ext
      lfilt2_ext = None
   #  
   if attfile == None: 
       raise IOError("The attitude file could not be found.") 
   return specfile, lfilt1, lfilt1_ext, lfilt2, lfilt2_ext, attfile     

def boresight(filter='uvw1',order=1,wave=260,
              r2d=77.0,date=0,chatter=0):
   ''' provide reference positions on the 
       UVOT filters for mapping and as function of 
       time for grisms. 
       
       This function name is for historical reasons, 
       and provides a key mapping function for the 
       spectral extraction.  
   
       The correct boresight of the (lenticular) filters 
       should be gotten from the Swift UVOT CALDB 
       as maintained by HEASARC. The positions here 
       are in some cases substantially different from
       the boresight in the CALDB. They are reference 
       positions for the spectral extraction algorithms 
       rather than boresight. 
       
       The grism boresight positions at 260nm (uv grism)
       and 420nm (visible grism) in first order are served
       in an uncommon format (in DET pixels) 
       by adding (77,77) to the lenticular filter 
       RAW coordinate.(see TELDEF file) the grism 
       boresight was measured in DET coordinates, 
       not RAW. (offset correction should be 104,78)

       Parameters
       ----------
       filter : str 
          one of {'ug200','uc160','vg1000','vc955',
          'wh','v','b','u','uvw1','uvm2','uvw2'}
       
       order : {0,1,2}
          order for which the anchor is needed

       wave : float
          anchor wavelength in nm

       r2d : float 
          additive factor in x,y to anchor position 

       date: long
          format in swift time (s)
          if 0 then provide the first order anchor 
          coordinates of the boresight for mapping 
          from the lenticular filter position 

       chatter : int 
          verbosity 

       Returns
       ------- 
       When *date* = 0:
       
       For translation: The boresight for a filter 
       (in DET pixels) by adding (77,77) to the 
       lenticular filter RAW coordinate (see TELDEF file)
       the grism boresight was measured in DET 
       (The default r2d=77 returns the correct 
       boresight for the grisms in detector 
       coordinates. To get the grism boresight in 
       detector image coordinates, subtract (104,78) 
       typically. The difference is due to the distortion
       correction from RAW to DET)
       
       When *date* is non-zero, and *order*=0:
       The zeroth order boresight  
      
          
       NOTE: 
       -----
       THE TRANSLATION OF LENTICULAR IMAGE TO GRISM 
       IMAGE IS ALWAYS THE SAME, INDEPENDENT OF THE 
       BORESIGHT.
       THEREFORE THE BORESIGHT DRIFT DOES NOT AFFECT 
       THE GRISM ANCHOR POSITIONS AS LONG AS THE DEFAULT 
       BORESIGHT POSITIONS ARE USED. 
       [Becase those were used for the calibration].

       However, the zeroth order "reference" position 
       drift affects the "uvotgraspcorr" - derived 
       WCS-S. The positions used 

       History: 
         2014-01-04 NPMK : rewrite to inter/extrapolate 
         the boresight positions
       
   '''
   from scipy.interpolate import interp1d
   import numpy as np
   
   filterlist = ['ug200','uc160','vg1000','vc955',
           'wh','v','b','u','uvw1','uvm2','uvw2']
   if filter == 'list': return filterlist
   grismfilters = ['ug200','uc160','vg1000','vc955']
   lenticular = ['v','b','u','uvw1','uvm2','uvw2']
   
   #old pixel offset anchor based on pre-2010 data
   # dates in swift time, drift [x.y] in pixels 
   #dates=[209952000,179971200,154483349,139968000,121838400]
   #drift=[ [0,0], [+2.4,-2.0], [+3.4,-3.0], [+6.4,-10], [+6.4,-10]]
   
   # data from Frank's plot (email 2 dec 2013, uvw1 filter)
   # original plot was in arcsec, but the drift converted 
   # to pixels. uvw1 seems representative (except for white)
   swtime = np.array([  
         1.25000000e+08,   1.39985684e+08,   1.60529672e+08,
         1.89248438e+08,   2.23489068e+08,   2.46907209e+08,
         2.66126366e+08,   2.79601770e+08,   2.89763794e+08,
         3.01251301e+08,   3.13180634e+08,   3.28423998e+08,
         3.43445470e+08,   3.59351249e+08,   3.75257678e+08,
         4.50000000e+08])
   boredx = (np.array([-1.6, -0.870,0.546,1.174,2.328,2.47,
        2.813,3.076,3.400,3.805,4.149,4.656,
        5.081,5.607,6.072,8.56 ])-1.9)/0.502
   boredy = (np.array([ -0.75,-2.197,-4.857,-6.527,
        -7.098,-7.252,-7.142,-7.560,
        -7.670,-8.000,-8.043,-8.395,
        -8.637,-9.142,-9.670,-11.9])+6.8)/0.502
   # I assume the same overall drift for the grism 
   # boresight (in pixels). Perhaps a scale factor for the 
   # grism would be closer to 0.56 pix/arcsec 
   # the range has been extrapolated for better interpolation
   # and also to support the near future. The early
   # time extrapolation is different from the nearly constant
   # boresight in the teldef but within about a pixel.
   # I think the extrapolation is more accurate.
   fx = interp1d(swtime,boredx,bounds_error=False,fill_value="extrapolate")
   fy = interp1d(swtime,boredy,bounds_error=False,fill_value="extrapolate")
   
   # reference anchor positions          
   reference0 = {'ug200': [1449.22, 707.7],
                 'uc160': [1494.9 , 605.8], #[1501.4 , 593.7], # ?[1494.9, 605.8],
                 'vg1000':[1506.8 , 664.3],
                 'vc955': [1542.5 , 556.4]} 
                  
   # DO NOT CHANGE THE FOLLOWING VALUES AS THE WAVECAL DEPENDS ON THEM !!!
   reference1 = {'ug200': [ 928.53,1002.69],
                 'uc160': [1025.1 , 945.3 ], 
                 'vg1000':[ 969.3 ,1021.3 ],
                 'vc955': [1063.7 , 952.6 ]}            
                          
   if (filter in grismfilters):
      if (date > 125000000) and (order == 0):
          anchor = reference0[filter]
          anchor[0] += r2d-fx(date)
          anchor[1] += r2d-fy(date)
          return anchor
      elif (date > 125000000) and (order == 1):   
          anchor = reference1[filter]
          anchor[0] += r2d-fx(date)
          anchor[1] += r2d-fy(date)
          return anchor
      elif order == 1:    
          anchor = reference1[filter]
          anchor[0] += r2d
          anchor[1] += r2d
          return anchor
      elif order == 0:  
          raise RuntimeError(
          "The zeroth order reference position needs a date")  
      else:
          return reference1[filter]       
                  
   elif (date > 125000000) and (filter in lenticular):
      ref_lent = {'v':[951.74,1049.89],
                  'b':[951.87,1049.67],
                  'u':[956.98,1047.84],
                  'uvw1':[951.20,1049.36],
                  'uvm2':[949.75,1049.30],
                  'uvw2':[951.11,1050.18]}
      anchor = ref_lent[filter]
      anchor[0] += r2d-fx(date)
      anchor[1] += r2d-fy(date)
      return anchor
      
   elif (date > 122000000) and (filter == 'wh'):
      print("approximate static white filter boresight")
      if date > 209952000:
         return 949.902+r2d, 1048.837+r2d         
      elif date > 179971200:
         return 953.315+r2d, 1048.014+r2d        
      elif date > 154483349:
         return 954.506+r2d, 1043.486+r2d
      elif date > 139968000:
         return 956.000+r2d, 1039.775+r2d
      elif date >  121838400:
         return 956.000+r2d, 1039.775+r2d      
      else: return filterlist

   else:
      # this is the version used initially *(changed 2 june 2009)
      # DO NOT CHANGE THESE VALUES AS THE WAVECAL DEPENDS ON THEM !!!
      if   filter == 'uvw1': return 954.61+r2d, 1044.66+r2d
      elif filter == 'wh'  : return 954.51+r2d, 1043.49+r2d
      elif filter == 'v'   : return 955.06+r2d, 1045.98+r2d 
      elif filter == 'b'   : return 955.28+r2d, 1045.08+r2d 
      elif filter == 'u'   : return 960.06+r2d, 1043.33+r2d
      elif filter == 'uvm2': return 953.23+r2d, 1044.90+r2d 
      elif filter == 'uvw2': return 953.23+r2d, 1044.90+r2d
      elif filter == 'w1'  : return 954.61+r2d, 1044.66+r2d
      elif filter == 'm2'  : return 953.23+r2d, 1044.90+r2d 
      elif filter == 'w2'  : return 953.23+r2d, 1044.90+r2d
      elif filter == 'ug200':       
          if order == 1:
             if wave == 260: return 928.53+r2d,1002.69+r2d
      elif filter == 'uc160':       
          if order == 1:
             if wave == 260: return 1025.1+27+r2d,945.3+r2d
      elif filter == 'vg1000': 
          #elif order == 1: return 948.4+r2d, 1025.9+r2d
          if order == 1: return 969.3+r2d, 1021.3+r2d
      elif filter == 'vc955':
          if order == 1: return 1063.7+r2d, 952.6+r2d
   raise IOError("valid filter values are 'wh','v',"\
        "'b','u','uvw1','uvm2','uvw2','ug200',"\
        "'uc160','vg1000','vc955'\n")    

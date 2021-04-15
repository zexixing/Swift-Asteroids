import numpy as np
def spec_curvature(wheelpos,anchor,order=1,):
   '''Find the coefficients of the polynomial for the curvature.   
   
   Parameters
   ----------
   wheelpos : int, {160,200,955,1000}
      grism filter position in filter wheel
   anchor : list, array
      anchor position in detector coordinates (pixels)
   order : int
      the desired spectral order  
         
   Returns
   -------
      Provides the polynomial coefficients for y(x). 
   
   Notes
   -----
   The curvature is defined with argument the pixel coordinate in the dispersion 
   direction with reference to the the anchor coordinates in det-img 
   coordinates. The polynomial returns the offset normal to the dispersion.
   
   - 2011-03-07 Paul Kuin, initial version 
   - 2011-08-02 fixed nominal coefficients order=1
   '''
   from scipy import interpolate
   from numpy import array
   xin = anchor[0] -104
   yin = anchor[1]  -78
   if ((wheelpos == 1000) ^ (wheelpos == 955)):
      # return y = 0 + 0.0*x coefficient
      return array([0.,0.])

   elif wheelpos == 160:

     if order == 1:

       tck_c1= [array([0.,0.,0.,0.,2048.,  2048.,  2048.,  2048.]), \
          array([0.,0.,0.,0.,  2048.,  2048.,  2048.,  2048.]), \
          array([ 0.1329227 , -0.28774943,  0.13672294, -0.18436127, -0.19086855,\
          0.23071908, -0.21803703,  0.11983982,  0.16678715, -0.2004285 ,\
          0.12813155, -0.13855324, -0.1356009 ,  0.11504641, -0.10732287,\
          0.03374111]),3,3]
          
       tck_c2 = [array([0.,0.,0.,0.,  2048.,  2048.,  2048.,  2048.]),\
          array([0.,0.,0.,0.,  2048.,  2048.,  2048.,  2048.]),\
          array([ -3.17463632e-04,   2.53197376e-04,  -3.44611897e-04,\
         4.81594388e-04,   2.63206764e-04,  -3.03314305e-04,\
         3.25032065e-04,  -2.97050826e-04,  -3.06358032e-04,\
         3.32952612e-04,  -2.79473410e-04,   3.95150704e-04,\
         2.56203495e-04,  -2.34524716e-04,   2.75320861e-04,\
        -6.64416547e-05]),3,3]
       
       tck_c3 = [array([ 0.,0.,0.,0.,2048.,  2048.,  2048.,  2048.]),\
          array([ 0.,0.,0.,0.,2048.,  2048.,  2048.,  2048.]),\
          array([ -4.14989592e-07,   5.09851884e-07,  -4.86551197e-07,\
          1.33727326e-07,   4.87557866e-07,  -5.51120320e-07,\
          5.76975007e-07,  -3.29793632e-07,  -3.42589204e-07,\
          3.00002959e-07,  -2.90718693e-07,   5.57782883e-08,\
          2.20540397e-07,  -1.62674045e-07,   8.70230076e-08,\
         -1.13489556e-07]),3,3]
                     
       #coef = array([interpolate.bisplev(xin,yin,tck_c3),interpolate.bisplev(xin,yin,tck_c2),\
       #              interpolate.bisplev(xin,yin,tck_c1), 0.])
       coef = array([interpolate.bisplev(xin,yin,tck_c3)*0.5,interpolate.bisplev(xin,yin,tck_c2)*0.5,\
                     interpolate.bisplev(xin,yin,tck_c1)*0.5, 0.])   #~FIXME:      
       return coef
       
     elif order == 2: 
        tck_c0 = [array([ 0., 0., 0., 0., 1134.78683, 2048., 2048., 2048., 2048.]), \
                  array([ 0., 0., 0., 0., 871.080060, 2048., 2048., 2048., 2048.]), \
        array([-110.94246902,   15.02796289,  -56.20252149,  -12.04954456,\
        311.31851187,  -31.09148174,  -48.44676102,   85.82835905,\
        -73.06964994,   99.58445164,   46.47352776,   11.29231744,\
        -68.32631894,   88.68570087,  -34.78582366,  -33.71033771,\
          6.89774103,   25.59082616,   23.37354026,   49.61868235,\
       -438.17511696,  -31.63936231,   28.8779241 ,   51.03055925,\
         16.46852299]), 3, 3]

        tck_c1 = [array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
                  array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([ 0.52932582, -0.76118033,  0.38401924, -0.189221  , -0.45446129,\
        0.73092481, -0.53433133,  0.12702548,  0.21033591, -0.45067611,\
        0.32032545, -0.25744487, -0.06022942,  0.22532666, -0.27174491,\
        0.03352306]), 3, 3]

        tck_c2 = [array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
                  array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([ -4.46331730e-04,   3.94044533e-04,  -1.77072490e-04,\
         2.09823843e-04,   3.02872440e-04,  -6.23869655e-04,\
         5.44400661e-04,  -3.70038727e-04,  -1.60398389e-04,\
         4.90085648e-04,  -4.91436626e-04,   4.62904236e-04,\
         4.05692472e-05,  -2.34521165e-04,   3.04866621e-04,\
        -1.25811263e-04]), 3, 3]

       #tck_c0 = [array([0.,0.,  1132.60995961,  2048.,2048.]),
       #          array([0.,0.,   814.28303687,  2048.,2048.]),
       #          array([-49.34868162,  -0.22692399, -11.06660953,   5.95510567,
       #            -3.13109456,  37.63588808, -38.7797533 ,  24.43177327,  43.27243297]),1,1]
       #tck_c1 = [array([    0.,     0.,  2048.,  2048.]), 
       #          array([    0.,     0.,  2048.,  2048.]),
       #          array([ 0.01418938, -0.06999955, -0.00446343, -0.06662488]),1,1]
       #tck_c2 = [array([    0.,     0.,  2048.,  2048.]),
       #          array([    0.,     0.,  2048.,  2048.]), 
       #         array([ -9.99564069e-05, 8.89513468e-05, 4.77910984e-05, 1.44368445e-05]),1,1]
       
        coef = array([interpolate.bisplev(xin,yin,tck_c2),interpolate.bisplev(xin,yin,tck_c1),\
                     interpolate.bisplev(xin,yin,tck_c0)])
        return coef
       
     elif order == 3: 
       # not a particularly good fit.
       tck_c0 =   [array([0.,     0.,  1101.24169141,  2048.,2048.]), 
                   array([0.,     0.,   952.39879838,  2048.,2048.]), 
                   array([ -74.75453915,    7.63095536, -131.36395787,   11.14709189,
                            -5.52089337,   73.59327202,  -57.25048374,   37.8898465 ,
                            65.90098406]), 1, 1]          
       tck_c1 = [array([    0.,     0.,  2048.,  2048.]), 
                 array([    0.,     0.,  2048.,  2048.]), 
                 array([-0.04768498, -0.02044308,  0.02984554, -0.04408517]), 1, 1]
 
       coef = array([interpolate.bisplev(xin,yin,tck_c1),interpolate.bisplev(xin,yin,tck_c0)])             
       return coef
       
     elif order == 0:
       tck_c0 =         [array([    0.,     0.,  1075.07521348,  2048. ,2048.]),
                  array([    0.,     0.,  1013.70915889,  2048. ,2048.]),
                  array([ 130.89087966,   25.49195385,    5.7585513 ,  -34.68684878,
                          -52.13229007, -168.75159696,  711.84382717, -364.9631271 ,
                          374.9961278 ]),1,1]
       tck_c1 =         [array([    0.,     0.,  2048.,  2048.]),
                  array([    0.,     0.,  2048.,  2048.]),
                  array([ 0.08258587, -0.06696916, -0.09968132, -0.31579981]),1,1]
                  
       coef = array([interpolate.bisplev(xin,yin,tck_c1),interpolate.bisplev(xin,yin,tck_c0)])            
       return  coef
     else: 
       raise (ValueError)    

   elif wheelpos == 200:
   
     if order == 1:
        tck_c1 = [array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([-0.00820665, -0.06820851,  0.04475057, -0.06496112,  0.062989  , \
        -0.05069771, -0.01397332,  0.03530437, -0.17563673,  0.12602437,\
        -0.10312421, -0.02404978,  0.06091811, -0.02879142, -0.06533121,\
         0.07355998]), 3, 3]
        
        tck_c2 = [array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([  1.69259046e-04,  -1.67036380e-04,  -9.95915869e-05, \
         2.87449321e-04,  -4.90398133e-04,   3.27190710e-04, \
         2.12389405e-04,  -3.55245720e-04,   7.41048332e-04, \
        -4.68649092e-04,  -1.11124841e-04,   6.72174552e-04, \
        -3.26167775e-04,   1.15602175e-04,   5.78187743e-04, \
        -8.79488201e-04]), 3, 3]

        tck_c3 = [array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([    0.,     0.,     0.,     0.,  2048.,  2048.,  2048.,  2048.]),\
        array([  1.11106098e-07,   2.72305072e-07,  -7.24832745e-07,\
         4.65025511e-07,  -2.35416547e-07,  -3.87761080e-07,\
         1.05955881e-06,  -6.46388216e-07,   3.15103869e-07,\
         5.48402086e-07,  -1.44488974e-06,   6.52867676e-07,\
         1.14004672e-08,  -9.48879026e-07,   1.64082320e-06,\
        -8.07897628e-07]), 3, 3]

        # the linear fit fails at the right side (57020002) but is quite good otherwise:
        #tck_c1 = [array([    0.,     0.,  2048.,  2048.]), array([    0.,     0.,  2048.,  2048.]),\
        #          array([-0.02212781, -0.00873168, -0.00377861, -0.02478484]), 1, 1]
        # 
        #tck_c2 = [array([    0.,     0.,  2048.,  2048.]), array([    0.,     0.,  2048.,  2048.]),\
        #          array([ -6.75189230e-05,   6.19498966e-05,   5.22322103e-05, 7.75736030e-05]), 1, 1]
        #
        #tck_c3 = [array([    0.,     0.,  2048.,  2048.]), array([    0.,     0.,  2048.,  2048.]), \
        #          array([ -1.75056810e-09,  -3.61606998e-08,  -6.00321832e-09, -1.39611943e-08]), 1, 1] 
        coef = array([interpolate.bisplev(xin,yin,tck_c3),interpolate.bisplev(xin,yin,tck_c2),\
                     interpolate.bisplev(xin,yin,tck_c1), 0.])
        return coef
       
     elif order == 2: 
     
       tck_c0 = [array([0.,0.,   956.25596245,  2048.,2048.]),
           array([0.,0.,  1067.40622524,  2048.,2048.]),
           array([ 17.82135471,  -4.93884392,  20.55439437, -18.22869669,
        13.11429182,  41.2680039 ,   9.8050793 ,  32.72362507,  -6.56524782]), 1, 1]
        
       tck_c1 =  [array([    0.,     0.,  2048.,  2048.]),
           array([    0.,     0.,  2048.,  2048.]),
           array([ 0.02362119, -0.03992572,  0.0177935 , -0.10163929]),1, 1]
           
       tck_c2 =  [array([    0.,     0.,  2048.,  2048.]),
           array([    0.,     0.,  2048.,  2048.]),
           array([ -6.32035759e-05,   5.28407967e-05,  -8.87338917e-06, 8.58873870e-05]),1,1]
       coef = array([interpolate.bisplev(xin,yin,tck_c2),interpolate.bisplev(xin,yin,tck_c1),\
                     interpolate.bisplev(xin,yin,tck_c0)])
       return coef
       
     elif order == 3:  
     
       tck_c0 = [array([    0.        ,     0.        ,   807.44415249,  2048.,2048.]),
                  array([    0.        ,     0.        ,  1189.77686531,  2048.,2048.]),
                  array([-5436.10353688,   218.93823252,  -254.71035527,   -24.35684969,
                   23.26131493,    51.66273635,    37.89898456,    46.77095978,
                   63.22039872]), 1, 1]

       tck_c1 = [array([    0.,     0.,  2048.,  2048.]),
                 array([    0.,     0.,  2048.,  2048.]),
                 array([-0.02591263, -0.03092398,  0.00352404, -0.01171369]), 1, 1]
       coef = array([interpolate.bisplev(xin,yin,tck_c1),interpolate.bisplev(xin,yin,tck_c0)])             
       return coef
       
     elif order == 0:
       tck_c0 = [array([0.,0.,   798.6983833,  2048.,  2048.]),
                 array([0.,0.,  1308.9171309,  2048.,  2048.]),
                 array([ 1244.05322027,    24.35223956,  -191.8634177 ,  -170.68236661,
          -4.57013926, 20.35393124, -365.28237355,  -235.44828185, -2455.96232688]), 1, 1]
       tck_c1 =  [array([    0.,     0.,  2048.,  2048.]),
                  array([    0.,     0.,  2048.,  2048.]),
                  array([ 0.54398146, -0.04547362, -0.63454342, -0.49417562]),1,1]

       coef = array([interpolate.bisplev(xin,yin,tck_c1),interpolate.bisplev(xin,yin,tck_c0)])            
       return  coef

     else: 
       raise (ValueError)    
       
   else:
      print('spec_curvature: illegal wheelpos value')
      raise (ValueError)   

def readFluxCalFile(wheelpos,anchor=None,option="default",spectralorder=1,
    arf=None,msg="",chatter=0):
   """Read the new flux calibration file, or return None.
   
   Parameters
   ----------
   wheelpos : int, required
      the position of the filterwheel 
   
   kwargs: dict
    - **anchor** : list, optional
      coordinate of the anchor
      
    - **option** : str
      option for output selection: 
        option=="default" + anchor==None: old flux calibration
        option=="default" + anchor : nearest flux calibration + model extrapolation
        option=="nearest" : return nearest flux calibration
        option=="model" : model 
        
    - **spectralorder** : int
        spectral order (1, or 2)
        
    - **arf**: path     
        fully qualified path to a selected response file
        
    - **msg**: str
        buffer message list (to add to)    

   Returns
   -------
   None if not (yet) supported
   option == 'model' returns the (astropy/pyfits) fits HDU (header+data) of the model 
   option == 'nearest'
      returns the fits HDU of the nearest calfile
   option == 'default' and anchor == None:
      returns the fits HDU of the nearest calfile 
   option == 'default' and anchor position given (in detector coordinates) 
      returns the fits HDU and an 
      interpolating function fnorm(wave in A) for the flux correction
   msg : string comments  separated by \n 
   
   Notes
   -----        
                 
   2013-05-05 NPMKuin
   """
   try:  
     from astropy.io import fits
   except:
     import pyfits as fits
   import os 
   import sys
   import numpy as np
   from scipy import interpolate

   typeNone = type(None)
   grismname = "UGRISM"
   if wheelpos > 500: grismname  = "VGRISM"
   if (type(anchor) != typeNone):
      if (len(anchor) != 2):
         sys.stderr.write("input parameter named anchor is not of length 2")
      elif type(anchor) == str: 
         anchor = np.array(anchor, dtype=float)  

   check_extension = False
   # here the "latest" version of the calibration files has been hardcoded    
   # latest update:   
   if spectralorder == 1: 
          if wheelpos == 200:          
             calfile = 'swugu0200_20041120v105.arf'
             extname = "SPECRESPUGRISM200"
             model   = "ZEMAXMODEL_200"
          elif wheelpos == 160:
             calfile = 'swugu0160_20041120v105.arf'
             extname = "SPECRESPUGRISM160"
             model   = "ZEMAXMODEL_160"
          elif wheelpos == 955: 
             calfile = 'swugv0955_20041120v104.arf'
             extname = "SPECRESPVGRISM0955"
             model   = "ZEMAXMODEL_955"
          elif wheelpos == 1000: 
             calfile = 'swugv1000_20041120v105.arf'
             extname = "SPECRESPVGRISM1000"
             model   = "ZEMAXMODEL_1000"
          else:   
             raise RuntimeError( "FATAL: [uvotio.readFluxCalFile] invalid filterwheel position encoded" )
             
   elif spectralorder == 2:          
          # HACK: force second order to 'nearest' option 2015-06-30 
          option == "nearest"
          check_extension = True
          if wheelpos == 200:          
             calfile =  'swugu0200_2_20041120v999.arf' #'swugu0200_20041120v105.arf'
             extname = "SPECRESP0160GRISM2NDORDER"
             model   = ""
          elif wheelpos == 160:
             calfile = 'swugu0160_2_20041120v999.arf'  #'swugu0160_20041120v105.arf'
             extname = "SPECRESP0160GRISM2NDORDER"
             model   = ""
          elif wheelpos == 955: 
             calfile = 'swugv0955_2_20041120v999.arf' #'swugv0955_20041120v104.arf'
             extname = "SPECRESPVGRISM955"
             model   = ""
          elif wheelpos == 1000: 
             calfile = 'swugv1000_2_20041120v999.arf'  #swugv1000_20041120v105.arf'
             extname = "SPECRESPVGRISM1000"
             model   = ""
          else:   
             raise RuntimeError( "FATAL: [uvotio.readFluxCalFile] invalid filterwheel position encoded" )
          if chatter > 3:
             print("[uvotio.readFluxCalFile] "+calfile)
             print("       extname="+extname)
             print("       model="+model+"|")  
   else:         
             raise RuntimeError("spectral order not 1 or 2 - no effective area available")
             
             
   if chatter > 1:
      print("uvotio.readFluxCalFile attempt to read effective area file: ")
   if arf != None:
      if arf.upper() == "CALDB":
   # try to get the file from the CALDB
         os.getenv("CALDB")
         command="quzcif swift uvota - "+grismname+\
          " SPECRESP now now  wheelpos.eq."+str(wheelpos)+" > quzcif.out"
         os.system(command) 
         f = open("quzcif.out")
         records = f.readlines()
         f.close()
         os.system("rm -f quzcif.out")
         arf, extens = records[0].split()  
         arf = CALDB + "/data/swift/uvota/cpf/arf/"+arf     
         hdu = fits.open(arf)
       
      else:
      # path to arf is supplied
      # the format must give the full path (starting with "/" plus FITS extension
      # if no extension was supplied and there is only one, assume that's it.
      # check version==2, using presence of CBD70001 keyword and see if spectral order is right
         if chatter > 3: print(arf)
         try:  # get extension from path 
            if len(arf.split("+") ) == 2: 
               file, extens = arf.split("+")
            elif len(arf.split("[") ) == 2:
               file = arf.split("[")[0]
               extens = arf.split("[")[1].split("]")[0] 
            else:
               check_extension = True
            arf = file
         except: 
            raise IOError("The supplied effective area file name "+arf+" cannot be understood.")                   
       
         hdu = fits.open(arf)
         if check_extension:  # old version file 
            if hdu[1].header['CBD60001'].split("(")[1].split(")")[0] != spectralorder: 
               raise IOError("The supplied effective area file is not correct spectral order.")
            if ("CBD70001" not in hdu[extens].header) :  # old version
               print("Spectral order = %i. \t"%(spectralorder))
               print("Using the oldest version of the effective area. \n"+\
                    "Flux, COI correction will be wrong.")
               return hdu[extname],msg
   
   else:    # argument arf = None      
       uvotpy = os.getenv("UVOTPY") + '/uvotpy'
       arf = os.path.join(uvotpy,"calfiles",calfile)
       
       #try:    
       hdu = fits.open(arf)
       if check_extension:  # old version file 
            if hdu[1].header['CBD60001'].split("(")[1].split(")")[0] != spectralorder: 
               #raise IOError("The supplied effective area file is not correct spectral order.")
               print("Spectral oder = %i:\t"%(spectralorder))
               print("The supplied effective area file %s is not \n   for the correct spectral order."%(arf))
            if ("CBD70001" not in hdu[extname].header) :  # old version
               print("Using the oldest version of the effective area. \n"+\
                    "Flux, COI correction will be wrong.")
               return hdu[extname],msg
       #except:
       #   print "UVOTPY environment variable not set or calfiles directory entries missing" 
       #   pass      
       #   return None, msg   
      
   if chatter > 0: print("Spectral order = %i: using flux calibration file: %s"%(spectralorder,arf))   
   if chatter > 2: hdu.info()
   msg += "Flux calibration file: %s\n"%(arf.split('/')[-1])
   
   if (option == "default") | (option == "nearest"):
      
      if type(anchor) == typeNone:  # assume centre of detector
         anchor = [1000,1000]
      else:
         if (option == "default"): modelhdu = hdu[model]
         if wheelpos < 500:
            n2 = 16
         else: n2 = 12   
         names = []      # names extensions
         calanchors = [] # anchor positions for each calibration extension
         dist = []       # distances 
         for h in range(1,len(hdu)):
            N = hdu[h].header["EXTNAME"].upper()
            NL = N.split("_")
            if (len(NL) == 3):
               if( int(NL[2][1]) == spectralorder): 
                  names.append(N)
                  root, ankxy, ord = NL
                  ankcal = ankxy.split("AY")
                  ankcal[0] = float(ankcal[0].split("AX")[1])
                  ankcal[1] = float(ankcal[1])
                  calanchors.append([ankcal[0],ankcal[1]])
                  dist.append( (ankcal[0]-anchor[0])**2+(ankcal[1]-anchor[1])**2 )
         # now we have a list of valid extnames, and anchor positions 
         # for the calibration file fits-extensions
         dist = np.array(dist) 
         k = np.where( dist == np.min(dist) )[0][0]      
         cal = hdu[names[k]]
         print("Nearest effective area is %s  - selected"%(names[k]))
         msg += "Selected nearest effective area FITS extension %s\n"%(names[k])
         if (option == "nearest"): 
            return cal, msg
         try:  
            if chatter > 4: 
               print("ReadFluxCalFile:      calanchor ", calanchors[k]) 
               print("ReadFluxCalFile:         anchor ", anchor)
               print("ReadFluxCalFile: MODEL  extname ", modelhdu.header['extname']) 
            modelcalflux = getZmxFlux (calanchors[k][0],calanchors[k][1],modelhdu)
            modelobsflux = getZmxFlux (anchor[0],anchor[1],modelhdu)
            q = np.isfinite(modelcalflux) & np.isfinite(modelobsflux) 
            w = 10.0*modelhdu.data['WAVE']
            if chatter > 4: 
              print("ReadFluxCalFile:         check:  ")
              print("ReadFluxCalFile:         w.shape ",w.shape)
              print("ReadFluxCalFile:            =784*",n2," ?")
            w = w.reshape(n2,784)[q,0]
            fn = modelobsflux[q]/modelcalflux[q]
            w1 = 1650.0
            f1 = 1.0      # was f1 = (fn[1]-fn[0])/(w[1]-w[0])*(w1-w[0]) + fn[0]
            n = len(w)+2
            x = np.zeros(n,dtype=float)
            y = np.zeros(n,dtype=float)
            x[0] = w1
            y[0] = f1
            x[1:-1] = w
            y[1:-1] = fn
            x[-1] = 7000.
            y[-1] = y[-2]
            y[ y < 0 ] = 0.0
            fnorm = interpolate.interp1d(x,y,bounds_error=False, fill_value=0.)     
            msg += "Flux corrected for variation over detector using model\n"
            return cal, fnorm, msg
         except RuntimeError:
             pass
             print("WARNING: Failure to use the model for inter/extrapolation of the calibrated locations.")
             print("         Using Nearest Eaafective Area File for the Flux calibration.")
             fnorm = interpolate.interp1d([1600,7000],[1.,1.],)
             return cal, fnorm, msg 
   elif option == "model":
       return hdu[model]
   else:
       raise RuntimeError( "invalid option passed to readFluxCalFile") 
'''
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from uvotpy.uvotio import getZmxFlux
x0,y0=(1000.99762536, 1135.91559724)

x = list(range(1,2048))
y = list(range(1,2048))
a = np.array([x,y])
#plt.imshow(spec_curvature(160,a,order=1,)[0],vmin=-5e-8,vmax=1e-10)
#plt.imshow(spec_curvature(160,a,order=1,)[2],vmin=-0.02,vmax=0.008)
plt.imshow(spec_curvature(160,a,order=1,)[1],vmin=-5e-5,vmax=5e-5)
plt.plot(x0, y0, 'kx',MarkerSize=4)
#plt.xlim(x0-15,x0+15)
#plt.ylim(y0-15,y0+15)
plt.title('coefficients 3')
plt.show()
'''

#x = list(range(975,1025))
#y = list(range(1110,1160))
#a = np.array([x,y])
#image = np.zeros((50,50))
#for i in x:
#    print(i)
#    for j in y:
#        image[i-975][j-1110] = readFluxCalFile(160,anchor=[j,i])[1]
'''
from scipy import interpolate
hdu, fnorm, msg =readFluxCalFile(160,anchor=[y0-25,x0+25])
w = list(0.5*(hdu.data['WAVE_MIN']+hdu.data['WAVE_MAX']) )
r = list(hdu.data['SPECRESP'])
w.reverse()
r.reverse()
specrespfunc = interpolate.interp1d( w, r, bounds_error=False, fill_value=np.NaN )
print(fnorm(3800))
print(specrespfunc(3800))
'''
#print(image)
#plt.imshow(image)
#plt.plot(x0, y0, 'kx',MarkerSize=4)
#plt.show()


def getCalData(Xphi, Yphi, wheelpos,date, chatter=3,mode='bilinear',
   kx=1,ky=1,s=0,calfile=None,caldir=None, msg=''):
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
   th  = data.field('SP1SLOPE').reshape(N1,N1)
   if wheelpos == 955:
      #  first  order dispersion
      c10 = data.field('DISP1_0').reshape(N1,N1)
      c11 = data.field('DISP1_1').reshape(N1,N1)
      c12 = data.field('DISP1_2').reshape(N1,N1)
      c13 = data.field('DISP1_3').reshape(N1,N1)
      c14 = np.zeros(N1*N1).reshape(N1,N1)
      c1n = data.field('DISP1_N').reshape(N1,N1)
      #  second order 
      xp2 = data.field('DETX2ANK').reshape(N1,N1)
      yp2 = data.field('DETY2ANK').reshape(N1,N1) 
      c20 = data.field('DISP2_0').reshape(N1,N1)
      c21 = data.field('DISP2_1').reshape(N1,N1)
      c22 = data.field('DISP2_2').reshape(N1,N1)
      c2n = data.field('DISP2_N').reshape(N1,N1)
   else:
      #  first  order dispersion
      c10 = data.field('disp1_0').reshape(N1,N1)
      c11 = data.field('disp1_1').reshape(N1,N1)
      c12 = data.field('disp1_2').reshape(N1,N1)
      c13 = data.field('disp1_3').reshape(N1,N1)
      c14 = data.field('disp1_4').reshape(N1,N1)
      c1n = data.field('disp1_N').reshape(N1,N1)
      #  second order 
      xp2 = data.field('detx2ank').reshape(N1,N1)
      yp2 = data.field('dety2ank').reshape(N1,N1) 
      c20 = data.field('disp2_0').reshape(N1,N1)
      c21 = data.field('disp2_1').reshape(N1,N1)
      c22 = data.field('disp2_2').reshape(N1,N1)
      c2n = data.field('disp2_n').reshape(N1,N1)
   #
   #  no transform here. but done to lookup array
   #
   rx, ry = Xphi, Yphi 
   #
   #  test if within ARRAY boundaries
   #
   xfp = xf[0,:]
   yfp = yf[:,0]
   if ((rx < min(xfp)) ^ (rx > max(xfp))):
      inXfp = False
   else:
      inXfp = True
   if ((ry < min(yfp)) ^ (ry > max(yfp))):
      inYfp = False
   else:
      inYfp = True         
   #
   #    lower corner (ix,iy)
   # 
   if inXfp :
      ix  = max( np.where( rx >= xf[0,:] )[0] ) 
      ix_ = min( np.where( rx <= xf[0,:] )[0] ) 
   else:
      if rx < min(xfp): 
         ix = ix_ = 0
         print("WARNING: point has xfield lower than calfile provides")
      if rx > max(xfp): 
         ix = ix_ = N1-1   
         print("WARNING: point has xfield higher than calfile provides")
   if inYfp :   
      iy  = max( np.where( ry >= yf[:,0] )[0] ) 
      iy_ = min( np.where( ry <= yf[:,0] )[0] ) 
   else:
      if ry < min(yfp): 
         iy = iy_ = 0
         print("WARNING: point has yfield lower than calfile provides")
      if ry > max(yfp): 
         iy = iy_ = 27   
         print("WARNING: point has yfield higher than calfile provides")
   if inYfp & inXfp & (chatter > 2): 
      print('getCalData.                             rx,         ry,     Xank,        Yank ')
      print(ix, ix_, iy, iy_)
      print('getCalData. gridpoint 1 position: ', xf[iy_,ix_], yf[iy_,ix_], xp1[iy_,ix_], yp1[iy_,ix_])
      print('getCalData. gridpoint 2 position: ', xf[iy ,ix_], yf[iy ,ix_], xp1[iy ,ix_], yp1[iy ,ix_])
      print('getCalData. gridpoint 3 position: ', xf[iy ,ix ], yf[iy ,ix ], xp1[iy ,ix ], yp1[iy ,ix ])
      print('getCalData. gridpoint 4 position: ', xf[iy_,ix ], yf[iy_,ix ], xp1[iy_,ix ], yp1[iy_,ix ])   
   #
   #  exception at outer grid edges: 
   #
   if ((ix == N1-1) ^ (iy == N1-1) ^ (ix_ == 0) ^ (iy_ == 0)):
           
     # select only coefficient with order 4 (or 3 for wheelpos=955)
     print("IMPORTANT:")
     print("\nanchor point is outside the calibration array: extrapolating all data") 

     try: 
      if wheelpos == 955 :
        # first order solution
        q4 = np.where( c1n.flatten() == 3 )
        xf = xf.flatten()[q4]
        yf = yf.flatten()[q4]
        xp1 = xp1.flatten()[q4]
        yp1 = yp1.flatten()[q4]
        th  = th.flatten()[q4]
        c10 = c10.flatten()[q4]
        c11 = c11.flatten()[q4]
        c12 = c12.flatten()[q4]
        c13 = c13.flatten()[q4]
        c14 = np.zeros(len(q4[0]))
        c1n = c1n.flatten()[q4]
        mode = 'bisplines'
        # second order solution only when at lower or right boundary
        if (ix == N1-1) ^ (iy == 0):
          q2 = np.where( c2n.flatten() == 2 )[0]
          xp2 = xp2.flatten()[q2]
          yp2 = yp2.flatten()[q2] 
          c20 = c20.flatten()[q2]
          c21 = c21.flatten()[q2]
          c22 = c22.flatten()[q2]
          c2n = c2n.flatten()[q2]
        else:
          N2 = N1/2
          xp2 = np.zeros(N2) 
          yp2 = np.zeros(N2) 
          c20 = np.zeros(N2)
          c21 = np.zeros(N2)
          c22 = np.zeros(N2)
          c2n = np.zeros(N2)
          
      else: 
        print('HHHHHIIII')
        q4 = np.where( c1n.flatten() == 4 )
        xf = xf.flatten()[q4]
        yf = yf.flatten()[q4]
        xp1 = xp1.flatten()[q4]
        yp1 = yp1.flatten()[q4]
        th  = th.flatten()[q4]
        c10 = c10.flatten()[q4]
        c11 = c11.flatten()[q4]
        c12 = c12.flatten()[q4]
        c13 = c13.flatten()[q4]
        c14 = np.zeros(len(q4[0]))
        c1n = c1n.flatten()[q4]
        xp2 = xp2.flatten()[q4]
        yp2 = yp2.flatten()[q4] 
        c20 = c20.flatten()[q4]
        c21 = c21.flatten()[q4]
        c22 = c22.flatten()[q4]
        c2n = c2n.flatten()[q4]
     
      # find the anchor positions by extrapolation
      anker  = np.zeros(2)
      anker2 = np.zeros(2)
      tck1x = interpolate.bisplrep(xf, yf, xp1, xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19,kx=3,ky=3,s=None) 
      tck1y = interpolate.bisplrep(xf, yf, yp1, xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19,kx=3,ky=3,s=None) 
      tck2x = interpolate.bisplrep(xf, yf, xp1, xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19,kx=3,ky=3,s=None) 
      tck2y = interpolate.bisplrep(xf, yf, yp1, xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19,kx=3,ky=3,s=None) 
     
      anker[0]  = xp1i = interpolate.bisplev(rx,ry, tck1x) 
      anker[1]  = yp1i = interpolate.bisplev(rx,ry, tck1y)      
      anker2[0] = xp2i = interpolate.bisplev(rx,ry, tck2x) 
      anker2[1] = yp2i = interpolate.bisplev(rx,ry, tck2y) 
      
      # find the angle  
      
      tck = interpolate.bisplrep(xf, yf, th,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      thi = interpolate.bisplev(rx,ry, tck)
      
      # find the dispersion
      
      tck = interpolate.bisplrep(xf, yf, c10,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      c10i = interpolate.bisplev(rx,ry, tck)
      tck = interpolate.bisplrep(xf, yf, c11,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      c11i = interpolate.bisplev(rx,ry, tck)
      tck = interpolate.bisplrep(xf, yf, c12,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      c12i = interpolate.bisplev(rx,ry, tck)
      tck = interpolate.bisplrep(xf, yf, c13,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      c13i = interpolate.bisplev(rx,ry, tck)
      tck = interpolate.bisplrep(xf, yf, c14,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
      c14i = interpolate.bisplev(rx,ry, tck)
      
      if ((ix == N1-1) ^ (iy == 0)):
         tck = interpolate.bisplrep(xf, yf, c20,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
         c20i = interpolate.bisplev(rx,ry, tck)
         tck = interpolate.bisplrep(xf, yf, c21,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
         c21i = interpolate.bisplev(rx,ry, tck)
         tck = interpolate.bisplrep(xf, yf, c22,xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=3,ky=3,s=None) 
         c22i = interpolate.bisplev(rx,ry, tck)
      else:
         c20i = c21i = c22i = np.NaN 
      if chatter > 2: 
            print('getCalData. bicubic extrapolation  ') 
            print('getCalData. first order anchor position = (%8.1f,%8.1f), angle theta = %7.1f ' % (xp1i,yp1i,thi ))
            print('getCalData. dispersion first  order = ',c10i,c11i,c12i,c13i,c14i)
            if c20i == NaN:
               print(" no second order extracted ")
            else:   
               print('getCalData. second order anchor position = (%8.1f,%8.1f) ' % (xp2i,yp2i))
               print('getCalData. dispersion second order = ', c20i,c21i, c22i)
     except:    
        print("failed - ABORTING")
        raise    
        return
   else: 
   #
   #  reduce arrays to section surrounding point
   #  get interpolated quantities and pass them on 
   # 
      if mode == 'bisplines':
      # compute the Bivariate-spline coefficients
      # kx = ky =  3 # cubic splines (smoothing) and =1 is linear
         task = 0 # find spline for given smoothing factor
      #  s = 0 # 0=spline goes through the given points
      # eps = 1.0e-6  (0 < eps < 1)
         m = N1*N1
         if chatter > 2: print('\n getCalData. splines ') 
         qx = qy = np.where( (np.isfinite(xrf.reshape(m))) & (np.isfinite(yrf.reshape(m)) ) )
         tck1 = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], xp1.reshape(m)[qx],xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s) 
         tck2 = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], yp1.reshape(m)[qx],xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s) 
         xp1i = interpolate.bisplev(rx,ry, tck1)
         yp1i = interpolate.bisplev(rx,ry, tck2)
         tck3 = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], th.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         thi  = interpolate.bisplev(rx,ry, tck3)
         xp2i = 0
         yp2i = 0
             
         if chatter > 2: print('getCalData. x,y,theta = ',xp1i,yp1i,thi, ' second order ', xp2i, yp2i)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c10.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c10i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c11.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c11i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c12.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c12i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c13.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c13i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c14.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c14i = interpolate.bisplev(rx,ry, tck)
         if chatter > 2: print('getCalData. dispersion first order = ',c10i,c11i,c12i,c13i,c14i)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c20.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c20i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c21.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c21i = interpolate.bisplev(rx,ry, tck)
         tck  = interpolate.bisplrep(xrf.reshape(m)[qx], yrf.reshape(m)[qy], c22.reshape(m),xb=-0.19,xe=+0.19,yb=-0.19,ye=0.19, kx=kx,ky=ky,s=s)
         c22i = interpolate.bisplev(rx,ry, tck)
         if chatter > 2: print('getCalData. dispersion second order = ', c20i,c21i, c22i)
      #
      if mode == 'bilinear':
         xp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), xp1 ,chatter=chatter)
         yp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), yp1 ,chatter=chatter)
         thi  = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), th  )# ,chatter=chatter)
         c10i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c10 )#,chatter=chatter)
         c11i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c11 )#,chatter=chatter)
         c12i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c12 )#,chatter=chatter)
         c13i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c13 )#,chatter=chatter)
         c14i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c14 )#,chatter=chatter)
         xp2i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), xp2 )#,chatter=chatter)
         yp2i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), yp2 )#,chatter=chatter)
         c20i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c20 )#,chatter=chatter)
         c21i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c21 )#,chatter=chatter)
         c22i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), c22 )#,chatter=chatter)
         print('HHHHHIIIIII')
         return xf,yf,xp1,yp1
         if chatter > 1: 
            print('getCalData. bilinear interpolation') 
            print('getCalData. first order anchor position = (%8.1f,%8.1f), angle theta = %7.1f ' % (xp1i,yp1i,thi ))
            print('getCalData. dispersion first  order = ',c10i,c11i,c12i,c13i,c14i)
            print('getCalData. second order anchor position = (%8.1f,%8.1f) ' % (xp2i,yp2i))
            print('getCalData. dispersion second order = ', c20i,c21i, c22i)
      if mode == 'interp2d':
         x1 = xf[0,:].squeeze()
         x2 = yf[:,0].squeeze()
         xp1i = interpolate.interp2d(x1,x2,xp1,kind='linear')
         #same as bisplines with s=0 and k=1
         return
                    
   C_1 = np.array([c14i,c13i,c12i,c11i,c10i])
   C_2 = np.array([c22i,c21i,c20i])
   # 
   # only theta for the first order is available 
   cal.close()
   anker  = np.array([xp1i,yp1i]) 
   anker2 = np.array([xp2i,yp2i]) 
   if chatter > 0: 
      print('getCalData. anker [DET-pix]   = ', anker)
      print('getCalData. anker [DET-img]   = ', anker - [77+27,77+1])
      print('getCalData. second order anker at = ', anker2, '  [DET-pix] ') 
   return anker, anker2, C_1, C_2, thi, data, msg


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
 
def phi_fit(t,rx,ry):
   xp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), xp1 ,chatter=chatter)
   yp1i = bilinear( rx, ry, xf[0,:].squeeze(), yf[:,0].squeeze(), yp1 ,chatter=chatter)
   return xp1i*(1-t)/2+yp1i*(1+t)/2

'''
Xphi = 0.0030912260029784156 
Yphi = -0.004113097476086955 
wheelpos = 160 
date1 = 324957563.89912
calfile = None 
chatter = 2
new_anker = [1134.4416060653666, 1008.4896966428131]
#anker, anker2, C_1, C_2, angle, calibdat, msg4 = getCalData(Xphi,Yphi,wheelpos, date1, \
#         calfile=calfile, chatter=chatter)
xf,yf,xp1,yp1 = getCalData(Xphi,Yphi,wheelpos, date1,calfile=calfile, chatter=chatter)
from scipy import optimize
coef, coef_err = optimize.curve_fit(phi_fit,[-1,1],new_anker,[Xphi,Yphi])
print(coef)
print(phi_fit(-1,0.00282986, -0.00292843),phi_fit(1,0.00282986, -0.00292843))
# [ 0.00282986 -0.00292843]
'''

def findInputAngle(RA,DEC,filestub, ext, wheelpos=200, 
       lfilter='uvw1', lfilter_ext=None, 
       lfilt2=None,    lfilt2_ext=None, 
       method=None, attfile=None, msg="",
       uvotgraspcorr_on=True,
       update_pnt=True,  
       catspec=None, indir='./', chatter=2):
   '''Find the angles along the X,Y axis for the target distance from the bore sight.
   
   Parameters
   ----------
   RA,DEC : float
      sky position, epoch J2000, decimal degrees
      
   filestub : str
      part of filename consisting of "sw"+`obsid`
      
   ext : int
      number of the extension
      
   kwargs : dict
   
      - **wheelpos** : int, {160,200,955,1000}      
        grism filter selected in filter wheel
         
      - **lfilter**, **lfilt2** : str, {'uvw2','uvm2','uvw1','u','b','v'}      
        lenticular filter name before and after grism exposure
        
      - **lfilter_ext**, **lfilt2_ext** : int   
        lenticular filter extension before and after grism exposure 
        
      - **method** : str, {'grism_only'}
        if set to `grism_only`, create a temporary header to compute the 
        target input angles, otherwise use the lenticular file image.
      
      - **attfile** : str, path 
        full path+filename of attitude file
        
      - **catspec** : path
        optional full path to catalog spec file to use with uvotgraspcorr
        
      - **indir** : str, path
        data directory path
        
      - **uvotgraspcorr_on** : bool
           enable/disable update of the WCS keywords in the grism file using uvotgraspcorr
           
      - **update_pnt** : bool      
           enable/disable and update to the WCS keywords in the grism file from the 
           attitude file, prior to running uvotgraspcorr (if enabled)
      - **chatter** : int
        verbosity

   Returns
   -------
   anker_as : array
      offset (DX,DY) in arcsec in DET coordinate system of the source 
      from the boresight
      needs to be converted to input rays by applying transform.
          
   anker_field : array
      offset(theta,phi) in degrees from the axis for 
      the input field coordinates for the zemax model lookup             
      
   tstart : float
      start time exposure (swift time in seconds)
      
   msg : string 
      messages   
   
   
   Notes
   -----
   Provided a combined observation is available in a lenticular filter and a grism 
   (i.e., they were aquired in the the same observation,) this routine determines the
   input angles from the boresight.  Assumed is that the grism and lenticular filter 
   image have the same extension. 
   
   If not a lenticular filter image was taken just before or after the grism exposure, 
   the input angles are determined from the grism aspect only. Before running this, 
   run uvotgrapcorr on the grism image when there is no lenticular filter to get a 
   better aspect solution.
         
   '''
   # 2017-05-17 an error was found in fits header read of the extension of a second filter
   #            which was introduced when converting to astropy wcs transformations
   # 2015-06-10 output the lenticular filter anchor position
   #            and fix deleted second lenticular filter 
   # 2015-07-16 changeover to astropy.wcs from ftools 
   # 2010-07-11 added code to move existing uvw1 raw and sky files out of the way and cleanup afterwards.
   # npkuin@gmail.com
   from astropy import wcs
   import numpy as np
   try:
      from astropy.io import fits
   except:
      import pyfits as fits

   from uvotwcs import makewcshdr 
   import os, sys
   
   __version__ = '1.2 NPMK 20170517 NPMK(MSSL)'

   msg = ""
   lenticular_anchors = {}
   
   if (chatter > 1):
      print("uvotgetspec.getSpec(",RA,DEC,filestub, ext, wheelpos, lfilter, lfilter_ext, \
       lfilt2,    lfilt2_ext, method, attfile, catspec, chatter,')')
   
   if ( (wheelpos == 160) ^ (wheelpos == 200) ): 
      gfile = indir+'/'+filestub+'ugu_dt.img'
   elif ( (wheelpos == 955) ^ (wheelpos == 1000) ): 
      gfile = indir+'/'+filestub+'ugv_dt.img'
   else: 
      sys.stderr.write("uvotgetspec.findInputAngle: \n\tThe wheelpos=%s is wrong! \n"+\
          "\tAborting... could not determine grism type\n\n"%(wheelpos)) 
      return   
      
   if ((lfilter == None) & (lfilt2 == None)) | (method == 'grism_only') : 
      lfilter = 'fk' 
      method == 'grism_only'
      lfilter_ext = 1
       
   uw1rawrenamed = False
   uw1skyrenamed = False
   if method == 'grism_only':
       if chatter > 1: print("grism only method. Creating fake lenticular uvw1 file for grism position")
       # test if there is already a uvw1 raw or sky file before proceeding
     
       if chatter > 2: 
           print('wheelpos ',wheelpos)
           print('attfile  ',attfile)
       wheelp1 = wheelpos
       rawfile = makewcshdr(filestub,ext,
                            attfile,
                            wheelpos=wheelp1,
                            indir=indir,
                            catspec=catspec,
                            uvotgraspcorr_on=uvotgraspcorr_on,
                            update_pnt=update_pnt,
                chatter=chatter) 
       # note that the path rawfile  = indir+'/'+filestub+'ufk_sk.img'
       tempnames.append(filestub)
       tempntags.append('fakefilestub')
   
   if lfilter_ext == None: 
       lfext = ext 
   else: 
       lfext = lfilter_ext  
   
   ffile = indir+'/'+filestub+'uw1_sk.img'
   if lfilter == 'wh'   : ffile = indir+'/'+filestub+'uwh_sk.img'
   if lfilter == 'u'    : ffile = indir+'/'+filestub+'uuu_sk.img'
   if lfilter == 'v'    : ffile = indir+'/'+filestub+'uvv_sk.img'
   if lfilter == 'b'    : ffile = indir+'/'+filestub+'ubb_sk.img'
   if lfilter == 'uvm2' : ffile = indir+'/'+filestub+'um2_sk.img'
   if lfilter == 'uvw2' : ffile = indir+'/'+filestub+'uw2_sk.img'
   if lfilter == 'fk'   : ffile = indir+'/'+filestub+'ufk_sk.img'

   hf = fits.getheader(ffile,lfext)   
   hg = fits.getheader(gfile,ext)

   # check for losses in grism image
   if (' BLOCLOSS' in hg):
       if float(hg['BLOCLOSS']) != 0: 
           print('#### BLOCLOSS = '+repr(hg['BLOCLOSS']))
           msg += "BLOCLOSS=%4.1f\n"%(hg['BLOCLOSS'])
   if ('STALLOSS' in hg):
       if (float(hg['STALLOSS']) != 0): 
           print('#### STALLOSS = '+repr(hg['STALLOSS']))
           msg += "STALLOSS=%4.1f\n"%(hg['STALLOSS'])
   if ('TOSSLOSS' in hg):
       if float(hg['TOSSLOSS']) != 0: 
           print('#### TOSSLOSS = '+repr(hg['TOSSLOSS']))
           msg += "TOSSLOSS=%4.1f\n"%(hg['TOSSLOSS'])
   tstart = hg['TSTART']

   if chatter > 1: print('grism exposure time = ',hg['EXPOSURE'],'  seconds')
   
   RA_PNT  = hg['RA_PNT']
   DEC_PNT = hg['DEC_PNT']
   PA_PNT  = hg['PA_PNT']   # roll angle
   time    = hg['TSTART']   # time observation
   ra_diff  = RA - RA_PNT
   dec_diff = DEC - DEC_PNT
   if ((ra_diff > 0.4) ^ (dec_diff > 0.4) ): 
       sys.stderr.write( 
       "\nWARNING: \n\tthe difference in the pointing from the header to the RA,DEC parameter is \n"+\
       "\tlarge delta-RA = %f deg, delta-Dec = %f deg\n\n"%(ra_diff,dec_diff))
       
   W1 = wcs.WCS(hf,)
   xpix_, ypix_ = W1.wcs_world2pix(RA,DEC,0)    
   W2 = wcs.WCS(hf,key='D',relax=True)    
   x1, y1 = W2.wcs_pix2world(xpix_,ypix_,0)
   print(xpix_, ypix_)
   print(RA,DEC)
   xpix_new,ypix_new=W2.wcs_world2pix(-0.4722625 ,  0.05284593,0)
   print(xpix_new,ypix_new) 
   print(W1.wcs_pix2world(xpix_new,ypix_new,0))
   
   RAs = repr(RA)
   DECs= repr(DEC)
   exts = repr(ext)
   lfexts = repr(lfext)

   # tbd - get random number for temp file name
   from os import getenv,system
   #system('echo '+RAs+'  '+DECs+' > radec.txt ' )

   CALDB = getenv('CALDB')
   if CALDB == '': 
       print('the CALDB environment variable has not been set')
       return None
   HEADAS = getenv('HEADAS')
   if HEADAS == '': 
       print('The HEADAS environment variable has not been set')
       print('That is needed for the uvot Ftools ')
       return None 
       
   #command = HEADAS+'/bin/uvotapplywcs infile=radec.txt outfile=skyfits.out wcsfile=\"'\
   #          +ffile+'['+lfexts+']\" operation=WORLD_TO_PIX chatter='+str(chatter)
   #if chatter > 0: print command
   #system( command )

   #f = open('skyfits.out', "r")
   #line = f.read()
   #if chatter > 1: print 'skyfits.out: '+line
   #x1, y1 = (line.split())[2:4]
   #f.close  
   #system( 'echo '+repr(x1)+'  '+repr(y1)+'  > skyfits.in' )
   ## 
   #command = HEADAS+'/bin/uvotapplywcs infile=skyfits.in outfile=detmm.txt wcsfile=\"'\
   #          +ffile+'['+lfexts+']\" operation=PIX_TO_WORLD to=D chatter='+str(chatter)
   #if chatter > 1: print command
   #system( command )
   #f = open('detmm.txt', "r")
   #line = f.read()
   #if chatter > 1: print 'detmm: '+line
   #x1, y1 = line.split()[2:4]
   #f.close
   #x1 = float(x1)
   #y1 = float(y1)
   if chatter > 1:
       print("\t The [det]coordinates in mm are (%8.4f,%8.4f) " % ( x1, y1))
   # convert anchor in DET coordinate mm to pixels and arcsec from boresight
   anker_uvw1det = np.array([x1,y1])/0.009075+np.array((1100.5,1100.5))
   msg += "LFILT1_ANCHOR= [%6.1f,%6.1f]\n"%(anker_uvw1det[0],anker_uvw1det[1])
   lenticular_anchors.update({"lfilt1":lfilter,"lfilt1_anker":anker_uvw1det})
   
   if (x1 < -14) | (x1 > 14) | (y1 < -14) | (y1 > 14) :
      # outside detector 
      print("\nERROR: source position is not on the detector! Aborting...",(x1,y1))
      raise IOError("\nERROR: source position is not on the detector! ")
   
   if lfilter == "fk" : 
      l2filter = "uvw1"
   else: l2filter = lfilter   
   if wheelpos != 160:
       anker_uvw1det_offset = anker_uvw1det - np.array( boresight(filter=l2filter))  # use fixed default value boresight 
   else: 
       anker_uvw1det_offset = anker_uvw1det - np.array( boresight(filter=l2filter,date=209952100) )
       print(anker_uvw1det_offset)
       print(anker_uvw1det)
       print(np.array( boresight(filter=l2filter,date=209952100) ))
   Xphi, Yphi = anker_uvw1det_offset*0.502    
   as2deg = 1./3600.    
   print(Xphi,Yphi)
   print(Xphi*as2deg,Yphi*as2deg)

   # second lenticular filter 
   
   if lfilt2 != None: 
      if lfilt2 == 'wh'   : f2ile = indir+'/'+filestub+'uwh_sk.img'
      if lfilt2 == 'u'    : f2ile = indir+'/'+filestub+'uuu_sk.img'
      if lfilt2 == 'v'    : f2ile = indir+'/'+filestub+'uvv_sk.img'
      if lfilt2 == 'b'    : f2ile = indir+'/'+filestub+'ubb_sk.img'
      if lfilt2 == 'uvw1' : f2ile = indir+'/'+filestub+'uw1_sk.img'
      if lfilt2 == 'uvm2' : f2ile = indir+'/'+filestub+'um2_sk.img'
      if lfilt2 == 'uvw2' : f2ile = indir+'/'+filestub+'uw2_sk.img'
      if lfilt2 == 'fk'   : f2ile = indir+'/'+filestub+'ufk_sk.img'
      if lfilt2_ext == None: 
          lf2ext = ext 
      else: 
          lf2ext = lfilt2_ext  
      if chatter > 4: print("getting fits header for %s + %i\n"%(f2ile,lf2ext))
      hf2 = fits.getheader(f2ile,lf2ext)   
      W1 = wcs.WCS(hf2,)
      xpix_, ypix_ = W1.wcs_world2pix(RA,DEC,0)    
      W2 = wcs.WCS(hf2,key='D',relax=True)    
      x2, y2 = W2.wcs_pix2world(xpix_,ypix_,0)
      print('ZEXI XING:')    
      print({'xpix_':xpix_,'ypix_':ypix_,'x2':x2,'y2':y2})
      #command = HEADAS+'/bin/uvotapplywcs infile=radec.txt outfile=skyfits.out wcsfile=\"'\
      #       +f2ile+'['+str(lf2ext)+']\" operation=WORLD_TO_PIX chatter='+str(chatter)
      #if chatter > 0: print command
      #system( command )

      #f = open('skyfits.out', "r")
      #line = f.read()
      #if chatter > 1: print 'skyfits.out: '+line
      #x2, y2 = (line.split())[2:4]
      #f.close  
      #system( 'echo '+repr(x2)+'  '+repr(y2)+'  > skyfits.in' )
      # 
      #command = HEADAS+'/bin/uvotapplywcs infile=skyfits.in outfile=detmm.txt wcsfile=\"'\
      #       +f2ile+'['+str(lf2ext)+']\" operation=PIX_TO_WORLD to=D chatter='+str(chatter)
      #if chatter > 1: print command
      #system( command )
      #f = open('detmm.txt', "r")
      #line = f.read()
      #if chatter > 1: print 'detmm: '+line
      #x2, y2 = line.split()[2:4]
      #f.close
      #x2 = float(x1)
      #y2 = float(y1)
      if chatter > 2: 
          print(" The [det]coordinates in mm are (%8.4f,%8.4f) " % ( x2, y2))
      # convert anchor in DET coordinate mm to pixels and arcsec from boresight
      anker_lf2det = np.array([x2,y2])/0.009075+np.array((1100.5,1100.5))
      msg += "LFILT2_ANCHOR= [%6.1f,%6.1f]\n"%(anker_lf2det[0],anker_lf2det[1])
      lenticular_anchors.update({'lfilt2':lfilt2,'lfilt2_anker':anker_lf2det})
   
      if (x2 < -14) | (x2 > 14) | (y2 < -14) | (y2 > 14) :
         # outside detector 
         print("/nERROR: source position is not on the detector! Aborting...")
         raise IOError("/nERROR: source position in second lenticular filter is not on the detector! ")

   # combine lenticular filter anchors, compute (mean) offset, convert in units of degrees
   if lfilt2 != None:
      anker_uvw1det = (anker_uvw1det+anker_lf2det)*0.5 
   if lfilter == "fk" : 
      l2filter = "uvw1"
   else: l2filter = lfilter   
   if wheelpos != 160:
       anker_uvw1det_offset = anker_uvw1det - np.array( boresight(filter=l2filter))  # use fixed default value boresight 
   else: 
       anker_uvw1det_offset = anker_uvw1det - np.array( boresight(filter=l2filter,date=209952100) )
   Xphi, Yphi = anker_uvw1det_offset*0.502    
   as2deg = 1./3600.    

   # cleanup
   # taken out since file is needed still:
   #   if method == 'grism_only':  os.system('rm '+filestub+'uw1_??.img ')
   if uw1rawrenamed:      os.system('mv '+uw1newraw+' '+uw1oldraw)
   if uw1skyrenamed:      os.system('mv '+uw1newsky+' '+uw1oldsky)
   
   crpix = crpix1,crpix2 = hg['crpix1'],hg['crpix2']  
   crpix = np.array(crpix)   # centre of image
   cent_ref_2img = np.array([1100.5,1100.5])-crpix  
 
   if chatter > 4:
       sys.stderr.write('findInputAngle. derived undistorted detector coord source in lenticular filter 1 = (%8.5f,%8.5f)  mm '%(x1,y1))
       if lfilt2 != None:
           sys.stderr.write('findInputAngle. derived undistorted detector coord source in lenticular filter 2 = (%8.5f,%8.5f)  mm '%(x2,y2))
   if chatter > 2:  
       print('findInputAngle. derived undistorted detector coord lenticular filter 1         =  ',anker_uvw1det)
       print('findInputAngle. derived undistorted physical image coord lenticular filter 1   =  ',anker_uvw1det-cent_ref_2img)
       if lfilt2 != None:
          print('findInputAngle. derived undistorted detector coord lenticular filter 2         =  ',anker_lf2det)
          print('findInputAngle. derived undistorted physical image coord lenticular filter 1   =  ',anker_lf2det -cent_ref_2img)
       print('findInputAngle. derived boresight offset lenticular filter ',lfilter,' (DET pix): ',anker_uvw1det_offset)
       print('findInputAngle. derived boresight offset: (', Xphi, Yphi,') in \"  = (',Xphi*as2deg, Yphi*as2deg,') degrees')
   # cleanup temp files:   
   #system('rm radec.txt skyfits.out  skyfits.in detmm.txt')
   print(Xphi,Yphi)
   return Xphi*as2deg, Yphi*as2deg, tstart, msg, lenticular_anchors

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

'''
RA=166.61502
DEC=8.7495
filestub = 'sw'+'00091026003'
ext = int(1)
attfile = '/Users/zexixing/Research/swiftASTER/data/juno/00091026003/uvot/image/../../auxil/sw00091026003pat.fits'
indir='/Users/zexixing/Research/swiftASTER/data/juno/00091026003/uvot/image'

Xphi, Yphi, date1, msg3, lenticular_anchors = findInputAngle( 166.61502, 8.7495, filestub, ext,
      uvotgraspcorr_on=True, update_pnt=True,  msg="", \
      wheelpos=160, lfilter='uvw2', lfilter_ext=int(1), \
      lfilt2=None, lfilt2_ext=None, method=None, \
      attfile=attfile, catspec=None, indir=indir, chatter=2)

print(Xphi,Yphi)
'''


def findBackground(extimg,background_lower=[None,None], background_upper=[None,None],yloc_spectrum=100, 
    smo1=None, smo2=None, chatter=2):
   '''Extract the background from the image slice containing the spectrum.
   
   Parameters
   ----------
   extimg : 2D array
      image containing spectrum. Dispersion approximately along x-axis.
   background_lower : list
      distance in pixels from `yloc_spectrum` of the limits of the lower background region.
   background_upper : list
      distance in pixels from `yloc_spectrum` of the limits of the upper background region.   
   yloc_spectrum : int
      pixel `Y` location of spectrum
   smo1 : float
      smoothing parameter passed to smoothing spline fitting routine. `None` for default.  
   smo2 : float
      smoothing parameter passed to smoothing spline fitting routine. `None` for default. 
   chatter : int
      verbosity
      
   Returns
   -------
   bg : float
      mean background 
   bg1, bg2 : 1D arrays
      bg1 = lower background; bg2 = upper background
      inherits size from extimg.shape x-xoordinate
   bgsig : float
      standard deviation of background  
   bgimg : 2D array
      image of the background constructed from bg1 and/or bg2   
   bg_limits_used : list, length 4
      limits used for the background in the following order: lower background, upper background    
   (bg1_good, bg1_dis, bg1_dis_good, bg2_good, bg2_dis, bg2_dis_good, bgimg_lin) : tuple
      various other background measures    
   
   Notes
   -----
   
   **Global parameter**
   
     - **background_method** : {'boxcar','splinefit'}

   The two background images can be computed 2 ways:
   
     1. 'splinefit': sigma clip image, then fit a smoothing spline to each 
         row, then average in y for each background region
     2. 'boxcar':   select the background from the smoothed image created 
         by method 1 below.
     3. 'sigmaclip': do sigma clipping on rows and columns to get column 
         profile background, then clip image and mask, interpolate over masked 
         bits.  
     
   extimg  is the image containing the spectrum in the 1-axis centered in 0-axis
   `ank` is the position of the anchor in the image 
      
   I create two background images:
    
         1. split the image strip into 40 portions in x, so that the background variation is small
            compute the mean 
            sigma clip (3 sigma) each area to to the local mean
            replace out-of-image pixels with mean of whole image (2-sigma clipped)
            smooth with a boxcar by the smoothing factor
         2. compute the background in two regions upper and lower 
            linearly interpolate in Y between the two regions to create a background image    
      
      bg1 = lower background; bg2 = upper background
      
      smo1, smo2 allow one to relax the smoothing factor in computing the smoothing spline fit
      
   History
   -------
   -  8 Nov 2011 NPM Kuin complete overhaul     
          things to do: get quality flagging of bad background points, edges perhaps done here?    
   -  13 Aug 2012: possible problem was seen of very bright sources not getting masked out properly 
          and causing an error in the background that extends over a large distance due to the smoothing.
          The cause is that the sources are more extended than can be handled by this method. 
          A solution would be to derive a global background     
   -  30 Sep 2014: background fails in visible grism e.g., 57977004+1 nearby bright spectrum 
          new method added (4x slower processing) to screen the image using sigma clipping      
      '''
   import sys   
   import numpy as np   
   try:
     from convolve import boxcar
   except:
     from stsci.convolve import boxcar
   from scipy import interpolate  
   import stsci.imagestats as imagestats  
   import matplotlib.pyplot as plt  
     
   # initialize parameters
   cval = -1.0123456789
   background_method = 'boxcar'
   bgimg    = extimg.copy()
   out    = np.where( (np.abs(bgimg-cval) <= 1e-6) )
   in_img = np.where( (np.abs(bgimg-cval) >  1e-6) & np.isfinite(bgimg) )
   nx = bgimg.shape[1]  # number of points in direction of dispersion
   ny = bgimg.shape[0]  # width of the image     
   
   # sigma screening of background taking advantage of the dispersion being 
   # basically along the x-axis 
   _PROFILE_BACKGROUND_ = False
   if _PROFILE_BACKGROUND_:
      bg, u_x, bg_sig = background_profile(bgimg, smo1=30, badval=cval)
      u_mask = np.zeros((ny,nx),dtype=bool)
      for i in range(ny):
          u_mask[i,(bgimg[i,:].flatten() < u_x) & 
                np.isfinite(bgimg[i,:].flatten())] = True
                
      bkg_sc = np.zeros((ny,nx),dtype=float)
      # the following leaves larger disps in the dispersion but less noise; 
      # tested but not implemented, as it is not as fast and the mean results 
      # are comparable: 
      #for i in range(ny):
      #    uf = interpolate.interp1d(np.where(u_mask[i,:])[0],bgimg[i,u_mask[i,:]],bounds_error=False,fill_value=cval)
      #    bkg_sc[i,:] = uf(np.arange(nx))
      #for i in range(nx):
      #    ucol = bkg_sc[:,i]
      #    if len(ucol[ucol != cval]) > 0:
      #        ucol[ucol == cval] = ucol[ucol != cval].mean()   
      for i in range(nx):
          ucol = bgimg[:,i]
          if len(ucol[u_mask[:,i]]) > 0: 
              ucol[np.where(u_mask[:,i] == False)[0] ] = ucol[u_mask[:,i]].mean()
          bkg_sc[:,i] = ucol                        
      if background_method == 'sigmaclip':
          return bkg_sc  
      else:
      # continue now with the with screened image
          bgimg = bkg_sc    
   
   kx0 = 0 ; kx1 = nx # default limits for valid lower background  
   kx2 = 0 ; kx3 = nx # default limits for valid upper background  
   ny4 = int(0.25*ny) # default width of each default background region 
   
   sig1 = 1 # unit for background offset, width
   bg_limits_used = [0,0,0,0]  # return values used 

   ## in the next section I replace the > 2.5 sigma peaks with the mean  
   ## after subdividing the image strip to allow for the 
   ## change in background level which can be > 2 over the 
   ## image. Off-image parts are set to image mean. 
   # this works most times in the absence of the sigma screening,but 
   # can lead to overestimates of the background.  
   # the call to the imagestats package is only done here, and should 
   # consider replacement. Its not critical for the program.
   #   
   xlist = np.linspace(0,bgimg.shape[1],80)
   xlist = np.asarray(xlist,dtype=int)
   imgstats = imagestats.ImageStats(bgimg[in_img[0],in_img[1]],nclip=3)   
   bg = imgstats.mean
   bgsig  = imgstats.stddev
   
   if chatter > 2:
      sys.stderr.write( 'background statistics: mean=%10.2f, sigma=%10.2f '%
         (imgstats.mean, imgstats.stddev))
      
   # create boolean image flagging good pixels
   img_good = np.ones(extimg.shape,dtype=bool)
   # flag area out of picture as bad
   img_good[out] = False

   # replace high values in image with estimate of mean  and flag them as not good   
   
   for i in range(78):
      # after the sigma screening this is a bit of overkill, leave in for now 
      sub_bg = boxcar(bgimg[:,xlist[i]:xlist[i+2]] , (5,5), mode='reflect', cval=cval)
      sub_bg_use = np.where( np.abs(sub_bg - cval) > 1.0e-5 ) # list of coordinates
      imgstats = None
      if sub_bg_use[0].size > 0: 
         imgstats = imagestats.ImageStats(sub_bg[sub_bg_use],nclip=3)
         # patch values in image (not out of image) with mean if outliers
         aval = 2.0*imgstats.stddev
         img_clip_ = (
            (np.abs(bgimg[:,xlist[i]:xlist[i+2]]-cval) < 1e-6) | 
            (np.abs(sub_bg - imgstats.mean) > aval) | 
            (sub_bg <= 0.) | np.isnan(sub_bg) )
         bgimg[:,xlist[i]:xlist[i+2]][img_clip_] = imgstats.mean  # patch image
         img_good[:,xlist[i]:xlist[i+2]][img_clip_] = False       # flag patches 
   #plt.imshow(extimg,vmin=10,vmax=50)
   #plt.show()
   #plt.imshow(img_good)
   #plt.show()
   # the next section selects the user-selected or default background for further processing
   if chatter > 1: 
      if background_method == 'boxcar': 
         sys.stderr.write( "BACKGROUND METHOD: %s;  background smoothing = %s\n"%
             (background_method,background_smoothing))
      else:
         sys.stderr.write( "BACKGROUND METHOD:%s\n"(background_method ))
      
   if not ((background_method == 'splinefit') | (background_method == 'boxcar') ):
      sys.stderr.write('background method missing; currently reads : %s\n'%(background_method))

   if background_method == 'boxcar':        
      # boxcar smooth in x,y using the global parameter background_smoothing
      bgimg = boxcar(bgimg,background_smoothing,mode='reflect',cval=cval)
   
   if background_lower[0] == None:
      bg1 = bgimg[0:ny4,:].copy()
      bg_limits_used[0]=0
      bg_limits_used[1]=ny4
      bg1_good = img_good[0:ny4,:] 
      kx0 = np.min(np.where(img_good[0,:]))+10  # assuming the spectrum is in the top two thirds of the detector
      kx1 = np.max(np.where(img_good[0,:]))-10
   else:
      # no curvature, no second order:  limits 
      bg1_1= np.max(np.array([yloc_spectrum - sig1*background_lower[0],20 ]))
      #bg1_0=  np.max(np.array([yloc_spectrum - sig1*(background_lower[0]+background_lower[1]),0]))
      bg1_0=  np.max(np.array([yloc_spectrum - sig1*(background_lower[1]),0]))
      bg1 = bgimg[int(bg1_0):int(bg1_1),:].copy() 
      bg_limits_used[0]=bg1_0
      bg_limits_used[1]=bg1_1
      bg1_good = img_good[int(bg1_0):int(bg1_1),:] 
      kx0 = np.min(np.where(img_good[int(bg1_0),:]))+10  # assuming the spectrum is in the top two thirds of the detector   
      kx1 = np.max(np.where(img_good[int(bg1_0),:]))-10  # corrected for edge effects
      
   #if ((kx2-kx0) < 20): 
   #   print 'not enough valid upper background points'   

   if background_upper[0] == None:
      bg2 = bgimg[-ny4:ny,:].copy()
      bg_limits_used[2]=ny-ny4
      bg_limits_used[3]=ny
      bg2_good = img_good[-ny4:ny,:]
      kx2 = np.min(np.where(img_good[ny-1,:]))+10  # assuming the spectrum is in the top two thirds of the detector
      kx3 = np.max(np.where(img_good[ny-1,:]))-10
   else:   
      bg2_0= np.min(np.array([yloc_spectrum + sig1*background_upper[0],180 ]))
      #bg2_1=  np.min(np.array([yloc_spectrum + sig1*(background_upper[0]+background_upper[1]),ny]))
      bg2_1=  np.min(np.array([yloc_spectrum + sig1*(background_upper[1]),ny]))
      bg2 = bgimg[int(bg2_0):int(bg2_1),:].copy()
      bg_limits_used[2]=bg2_0
      bg_limits_used[3]=bg2_1
      bg2_good = img_good[int(bg2_0):int(bg2_1),:]
      kx2 = np.min(np.where(img_good[int(bg2_1),:]))+10  # assuming the spectrum is in the top two thirds of the detector
      kx3 = np.max(np.where(img_good[int(bg2_1),:]))-10
      
   #if ((kx3-kx2) < 20): 
   #   print 'not enough valid upper background points'   
      

   if background_method == 'boxcar': 
      bg1 = bg1_dis = bg1.mean(0)
      bg2 = bg2_dis = bg2.mean(0)
      bg1_dis_good = np.zeros(nx,dtype=bool)
      bg2_dis_good = np.zeros(nx,dtype=bool)
      for i in range(nx):
        bg1_dis_good[i] = np.where(bool(int(bg1_good[:,i].mean(0))))
        bg2_dis_good[i] = np.where(bool(int(bg2_good[:,i].mean(0))))
      
   if background_method == 'splinefit':  
   
      #  mean bg1_dis, bg2_dis across dispersion 
      
      bg1_dis = np.zeros(nx) ; bg2_dis = np.zeros(nx)
      for i in range(nx):
         bg1_dis[i] = bg1[:,i][bg1_good[:,i]].mean()
         if not bool(int(bg1_good[:,i].mean())): 
            bg1_dis[i] = cval      
         bg2_dis[i] = bg2[:,i][bg2_good[:,i]].mean()  
         if not bool(int(bg2_good[:,i].mean())): 
            bg2_dis[i] = cval
      
      # some parts of the background may have been masked out completely, so 
      # find the good points and the bad points   
      bg1_dis_good = np.where( np.isfinite(bg1_dis) & (np.abs(bg1_dis - cval) > 1.e-7) )
      bg2_dis_good = np.where( np.isfinite(bg2_dis) & (np.abs(bg2_dis - cval) > 1.e-7) )
      bg1_dis_bad = np.where( ~(np.isfinite(bg1_dis) & (np.abs(bg1_dis - cval) > 1.e-7)) )   
      bg2_dis_bad = np.where( ~(np.isfinite(bg2_dis) & (np.abs(bg2_dis - cval) > 1.e-7)) )
                 
      # fit a smoothing spline to each background 
         
      x = bg1_dis_good[0]
      s = len(x) - np.sqrt(2.*len(x)) 
      if smo1 != None: s = smo1
      if len(x) > 40: x = x[7:len(x)-7]  # clip end of spectrum where there is downturn
      w = np.ones(len(x))   
      tck1 = interpolate.splrep(x,bg1_dis[x],w=w,xb=bg1_dis_good[0][0],xe=bg1_dis_good[0][-1],k=3,s=s) 
      bg1 = np.ones(nx) *  (bg1_dis[x]).mean()  
      bg1[np.arange(kx0,kx1)] = interpolate.splev(np.arange(kx0,kx1), tck1)
   
      x = bg2_dis_good[0]
      s = len(x) - np.sqrt(2.*len(x))      
      if smo2 != None: s = smo1
      if len(x) > 40: x = x[10:len(x)-10]  # clip
      w = np.ones(len(x))   
      tck2 = interpolate.splrep(x,bg2_dis[x],w=w,xb=bg2_dis_good[0][0],xe=bg2_dis_good[0][-1],k=3,s=s)    
      bg2 = np.ones(nx) *  (bg2_dis[x]).mean()  
      bg2[np.arange(kx2,kx3)] = interpolate.splev(np.arange(kx2,kx3), tck2)
      
      # force bg >= 0:
      # spline can do weird things ?
      negvals = bg1 < 0.0
      if negvals.any(): 
         bg1[negvals] = 0.0
         if chatter > 1:
            print("background 1 set to zero in ",len(np.where(negvals)[0])," points")
      
      negvals = bg2 < 0.0
      if negvals.any(): 
         bg2[negvals] = 0.0
         if chatter > 1:
            print("background 1 set to zero in ",len(np.where(negvals)[0])," points")
   
   # image constructed from linear inter/extra-polation of bg1 and bg2
   
   bgimg_lin = np.zeros(ny*nx).reshape(ny,nx)
   dbgdy = (bg2-bg1)/(ny-1)
   for i in range(ny):
      bgimg_lin[i,:] = bg1 + dbgdy*i
      
   # interpolate background and generate smooth interpolation image 
   if ( (background_lower[0] == None) & (background_upper[0] == None)):
        # default background region
        dbgdy = (bg2-bg1)/150.0 # assuming height spectrum 200 and width extraction regions 30 pix each
        for i9 in range(bgimg.shape[0]):
           bgimg[i9,kx0:kx1] = bg1[kx0:kx1] + dbgdy[kx0:kx1]*(i9-25)
           bgimg[i9,0:kx0] = bg2[0:kx0]
           bgimg[i9,kx1:nx] = bg2[kx1:nx]
        if chatter > 2: print("1..BACKGROUND DEFAULT from BG1 and BG2")   
   elif ((background_lower[0] != None) & (background_upper[0] == None)):
     # set background to lower background region   
        for i9 in range(bgimg.shape[0]):
           bgimg[i9,:] = bg1 
        if chatter > 2: print("2..BACKGROUND from lower BG1 only")   
   elif ((background_upper[0] != None) & (background_lower[0] == None)):
     # set background to that of upper background region   
        for i9 in range(bgimg.shape[0]):
           bgimg[i9,:] = bg2
        if chatter > 2: print("3..BACKGROUND from upper BG2 only")   
   else:
     # linear interpolation of the two background regions  
        dbgdy = (bg2-bg1)/(background_upper[0]+0.5*background_upper[1]+background_lower[0]+0.5*background_lower[1]) 
        for i9 in range(bgimg.shape[0]):
           bgimg[i9,kx0:kx1] = bg1[kx0:kx1] + dbgdy[kx0:kx1]*(i9-int(100-(background_lower[0]+0.5*background_lower[1])))
           bgimg[i9,0:kx0] =  bg2[0:kx0]    # assuming that the spectrum in not in the lower left corner 
           bgimg[i9,kx1:nx] = bg2[kx1:nx]
        if chatter > 2: print("4..BACKGROUND from BG1 and BG2")   
      
   return bg, bg1, bg2, bgsig, bgimg, bg_limits_used, (bg1_good, bg1_dis, 
          bg1_dis_good, bg2_good, bg2_dis, bg2_dis_good, bgimg_lin)


import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from pathlib import Path
from astropy.nddata import CCDData
from astropy.visualization import hist
from convenience_functions import show_image
from astropy import units as u
import ccdproc as ccdp
from astropy.stats import sigma_clip
from astropy.stats.funcs import mad_std
from ccdproc import ImageFileCollection
import warnings
from astropy.utils.exceptions import AstropyWarning
from sys import exit
from ccdproc import ccd_process
from astropy.nddata import Cutout2D
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.optimize import curve_fit

print('HAVE FUN AND ENJOY!!!')	

print("DONE BY FELIPE.")

''' To get more in touch with the code and with processing images from CCD cameras, please look at https://docs.astropy.org/en/stable/index.html lo learn more about it. You can
    also look at the following tutorial for basic analysis which is really useful https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/00-00-Preface.html            '''



# Last modified: 11/05/2023





def get_master_bias(path, folder = False, compare = False, bins = 500, distribution = False, title = 'Bias Level Camera 04' ,warning = 'ignore'):
    e focus on detecti
    
    # This function generate a master bias frame by combining dark frames of the same exposure time.
    
    ''' path: path to the folder where the bias frames are located, make sure they have the same exposure time
        Example:  path = '/home/felipe/Hiwi/Pycode/folder_bias'
        THE CODE WONT WORK WITH FILE.TIF FORMAT, WORK WITH FILE.FIT FORMAT.                                '''
    
    ''' folder = if you set it True it will show what is inside your folder, by default it is False.       '''
    
    ''' compare = if you set it True it will plot a singlet bias and the combined bias to compare them.
                  By default it is False.                                                                  '''
    ''' warning = decide you want to ignore warning or not. By default it is ignore. If you want to
                  see them, set warning = 'default'                                                        '''
    ''' so far, for this code, the combine files cannot be already in the folder, otherwise errors will
        arise. For that, if there is a file with the HEADER = 'combined', the code will stop.  
        Make sure, a combined dark file is not in the folder.                                              '''
   
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    path1 = path
    path = Path(path) 
    
    im_collection = ImageFileCollection(path)   # read all the files in the folder given by path
    
    ##########################################################################################
    
    for im in im_collection.hdus(combined=True):  # seek if the combined bias is already there,
        if im != None:
            #raise Exception('COMBINED BIAS IS ALREADY IN THE FOLDER!')
            warnings.warn("COMBINED BIAS ALREADY IN THE FOLDER, IT WILL BE OVERWRITTEN!")
            #exit()                                 # if it is already there, quit the code
         
    #########################################################################################
    
    if folder == True:          # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(im_collection.summary)   # display the files in the notebook
    
    ##########################################################################################
    
    
    calibrated_biases = im_collection.files_filtered( object ='bias', include_path=True, combined =None )

    masterbias = ccdp.combine(calibrated_biases,
                             method='median',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6, unit='adu'
                            )

    masterbias.meta['combined'] = True
    
    bias_file_name = 'masterbias.fit'

    masterbias.write(path / bias_file_name, overwrite = True)
    
    
    #########################################################################################
    
    if compare == True:  # if true, will plot a single bias and the combined bias to compare them
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        show_image(CCDData.read(calibrated_biases[0],unit='adu').data, cmap='gray', ax=ax1, fig=fig)
        ax1.set_title('Single calibrated bias')
        show_image(masterbias.data, cmap='gray', ax=ax2, fig=fig)
        ax2.set_title('{} bias images combined'.format(len(calibrated_biases)))
        
    ##########################################################################################
    
    if distribution == True: # compare the distribution of 1 bias and the combined one
        
        print ('PLOTTING BIN DATA COULD TAKE SOME TIME.')
        
        im1 = CCDData.read(calibrated_biases[0],unit='adu').data.flatten()
        im2 = masterbias.data.flatten()   # we read the data using CCDData and flatten() to be 1-Dimensional
        
        m = np.mean(im1)  # for setting the limit of the plot
        s = np.std(im1)
        
        lower = m - 3*s
        upper = m + 3*s
        
        plt.figure(figsize=(9, 9))
        hist(im1, bins = bins, label='One single bias', alpha=1, histtype='step', range = (250, 600)) # histogram of the data
        hist(im2, bins = bins, label='{} Bias frames combined'.format(len(calibrated_biases)), alpha=1, histtype='step', range = (250, 600))
        plt.legend(fontsize = 13)
        plt.grid()
        plt.xlim(lower,upper)
        plt.xlabel('Count level in image',fontsize = 20)
        plt.ylabel('Number of pixels with that count', fontsize = 20)
        plt.title(title, fontsize = 20)
        plt.savefig( path1 + 'Hist_Bias_Counts.png', dpi = 500,  bbox_inches='tight');








def get_calibrated_dark(path_dark, path_master_bias , folder = False, compare = False, distribution = False , cmap = 'gray' , warning = 'ignore' ):
    
    # This function will calibrated our dark frames by subtracting bias frames.
    # we cannot just combined the dark frames, we need to calibrate them first, just as the scientific images.
    
    '''
       path_dark: path where your uncalibrated dark frames are.
       
       folder: if you set it True it will show what is inside your folder, by default it is False.
       
       warning: decide you want to ignore warning or not. By default it is ignore. If you want to
                 see them, set warning = 'default'.
       
       path_master_bias: path to the master bias obtained by combining bias frames.
       
       compare: if True, you will compare the dark before and after calibration. By default it is False.
        
       cmap: you can choose the colormap of the plot.
                  
       '''
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    #######################################################################################
    
    # here we modify the directory where the flat frames are so we create a new folder for the calibrated flats
    
    
    import shutil
    
    folder_path = Path(path_dark + 'calibrated_darks/')

    # Delete the existing folder and all its contents
    if folder_path.is_dir():
    	shutil.rmtree(folder_path)
    
    if path_dark[-1] == '/':                
        
        path_dark_new = path_dark + 'calibrated_darks/'
        path_calibrated_darks = Path(path_dark_new)
        path_calibrated_darks.mkdir()
        print('THE CALIBRATED DARKS ARE PLACE IN:   ' + path_dark_new)
    
    if path_dark[-1] != '/':
        
        path_dark_new = path_dark + '/calibrated_darks/'
        path_calibrated_darks = Path(path_dark_new)
        path_calibrated_darks.mkdir()
        print('THE CALIBRATED DARKS ARE PLACE IN:   ' + path_dark_new)
       
    
    #######################################################################################    
   
    path_dark1 = path_dark
    
    path_dark = Path(path_dark)
    
    darks = ImageFileCollection(path_dark)
    
    path_master_bias = Path(path_master_bias)
    
    master_bias = CCDData.read(path_master_bias)
   
    #######################################################################################

    if folder == True:
        
        # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(darks.summary)   # display the files in the notebook
        
    ########################################################################################
    
    #flat_times = set(flats.summary['exptime']) # the exposure times of your dark frames
    
    #flat_times_list = list(flat_times)
    
    #if len(flat_times_list) > 1:  # we are checking if there are different exposure times which should not be the case
        
    #    raise Exception('YOU HAVE FLAT WITH DIFFERENT EXPOSURE TIME, PLEASE CHECK IT.')
    #    print(f'EXPOSURE TIME OF YOUR FLAT FRAMES ARE {flat_times_list} SECONDS.')
    #    #exit()  
        
    #if len(flat_times_list) == 1:
    
    #    print(f'EXPOSURE TIME OF YOUR FLAT FRAMES IS {flat_times_list[0]} SECONDS.')
        
    #########################################################################################
    	
    
    #exp_time_master_dark = master_dark.header['exptime']
    
    #if exp_time_master_dark != flat_times_list[0]: # we are checking if flat and dark have the same exposure time.
        
    #    raise Exception('MASTER DARK AND FLATS FRAMES HAVE DIFFERENT EXPOSURE TIME, PLEASE CHECK IT. THEY SHOULD HAVE THE SAME EXPOSURE TIME.')
    #    #exit()
        
    #########################################################################################
    
    for ccd, file_name in darks.ccds( object = 'DARK' , combined = None,     
                                      ccd_kwargs={'unit': 'adu'} ,   # CCDData requires a unit for the image if 
                                      return_fname = True ):         # Provide the file name too.                         
                                                                                        
        # Subtract the dark current 
        
        ccd = ccdp.subtract_bias(ccd, master_bias)

        # Save the result; there are some duplicate file names so pre-pend "calibrated_"
        
        ccd.write( path_calibrated_darks / ('calibrated_' + file_name))
        
    ##########################################################################################
    
    # here we plot the calibrated and ucalibrated flat. I have chosen it to be the first one because none of them should be more special than the other ones.
    
    if compare == True:
            
            dark_calib = ImageFileCollection(path_calibrated_darks)
            
            im2 = dark_calib.files_filtered(object = 'dark', include_path = True)
            
            im1 = darks.files_filtered(object = 'dark' ,include_path=True)
            
            no_calibrated = CCDData.read(im1[0],unit='adu')
            
            calibrated = CCDData.read(im2[0],unit='adu')
            
            fig, axes = plt.subplots(1, 2, figsize=(20, 10))

            show_image(no_calibrated, cmap=cmap, fig=fig, ax=axes[0])
            axes[0].set_title('Uncalibrated Dark', fontsize = 20)

            show_image(calibrated, cmap= cmap, fig=fig, ax=axes[1])
            axes[1].set_title('Calibrated Dark', fontsize = 20)
            
    ###################################################################################################################

    if distribution == True: # compare the distribution of 1 dark and the combined one
        
            print ('PLOTTING BIN DATA COULD TAKE SOME TIME.')
            
            dark_calib = ImageFileCollection(path_calibrated_darks)
            
            im2 = dark_calib.files_filtered(object = 'dark', include_path = True)
            
            im1 = darks.files_filtered(object	 = 'dark' ,include_path=True)
            
            no_calibrated = CCDData.read(im1[0],unit='adu').data.flatten()
            
            calibrated = CCDData.read(im2[0],unit='adu').data.flatten()
        
            m1 = np.mean(no_calibrated)  # for setting the limit of the plot
            s1 = np.std(no_calibrated)
        
            lower1 = m1 - 1*s1
            upper1 = m1 + 1*s1
            
            m2 = np.mean(calibrated)  # for setting the limit of the plot
            s2 = np.std(calibrated)
        
            lower2 = m2 - 1*s2
            upper2 = m2 + 1*s2
            
            fig, axes = plt.subplots(1, 2, figsize=(20, 10))
            
            hist(no_calibrated, bins= 5500, label='Uncalibrated dark', ax=axes[0]) # histogram of the data
            axes[0].set_title('Uncalibrated Dark', fontsize = 20)
            axes[0].set_xlabel('Count level in image',fontsize = 20)
            axes[0].set_ylabel('Number of pixels with that count',fontsize = 20)
            axes[0].set_xlim([lower1, upper1])
            axes[0].grid()


            hist(calibrated, bins= 5500, label='Calibrated dark', ax=axes[1])
            axes[1].set_title('Calibrated Dark', fontsize = 20)
            axes[1].set_xlabel('Count level in image',fontsize = 20)
            axes[1].set_xlim([lower2, upper2])
            axes[1].grid()
            
            plt.savefig( path_dark1 + 'Hist_dark_Calibrated.png', dpi = 500,  bbox_inches='tight');
            
            #plt.xlim(lower,upper)
            #plt.xlabel('Count level in image',fontsize = 20)
            #plt.ylabel('Number of pixels with that count', fontsize = 20)
    









def get_master_dark(path, folder = False, compare = False, distribution = False, title = 'Master Dark Camera 04' , same_time = True ,warning = 'ignore'):
    
    
    # This function generate a master dark frame by combining dark frames of the same exposure time.
    
    ''' path: path to the folder where the dark frames are located, make sure they have the same exposure time
        Example:  path = '/home/felipe/Hiwi/Pycode/folder_dark'
        THE CODE WONT WORK WITH FILE.TIF FORMAT, WORK WITH FILE.FIT FORMAT.                                '''
    
    ''' folder = if you set it True it will show what is inside your folder, by default it is False.       '''
    
    ''' compare = if you set it True it will plot a singlet dark and the combined dark to compare them.
                  By default it is False. '''
                  
    ''' same_time: If True, it will combined only dark frames with same exposure time. If False it will combine regardless
                   the exposure time but they have to be calibrated before. By defaukt True. '''
                                                                                   
    ''' warning = decide you want to ignore warning or not. By default it is ignore. If you want to
                  see them, set warning = 'default'      
                                                                    '''
    ''' so far, for this code, the combine files cannot be already in the folder, otherwise errors will
        arise. For that, if there is a file with the HEADER = 'combined', the code will stop.  
        Make sure, a combined dark file is not in the folder.                                              '''
   
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    path1 = path
    path = Path(path) 
    
    im_collection = ImageFileCollection(path)   # read all the files in the folder given by path
    
    ##########################################################################################
    
    for im in im_collection.hdus(combined=True):  # seek if the combined dark is already there,
        if im != None:
            #raise Exception('COMBINED DARK IS ALREADY IN THE FOLDER!')
            warnings.warn("COMBINED DARK ALREADY IN THE FOLDER, IT WILL BE OVERWRITTEN!")
            #exit()                                 # if it is already there, quit the code
         
    #########################################################################################
    
    if folder == True:          # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(im_collection.summary)   # display the files in the notebook
    
    dark_times = list(set(im_collection.summary['exptime'])) # the exposure times of your dark frames
    
    print('EXPOSURE TIME OF YOUR DARK FRAMES ARE ' + str(dark_times)+ ' SECONDS.')
    
    ##########################################################################################
    
    if same_time == True:
    
        for exp_time in sorted(dark_times):  # this will combine all the dark frames into one
            calibrated_darks = im_collection.files_filtered(object = 'dark', exptime=exp_time,include_path=True, combined = None)
            masterdark = ccdp.combine(calibrated_darks, 
                                         method='median',
                                         sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                         sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                                         mem_limit=350e6)
        # ccdp.combine combine the darks                      
        masterdark.meta['combined'] = True  # add to the combined file the HEADER combined = 'True'
 
        dark_file_name = 'masterdark{:6.3f}.fit'.format(exp_time)  # set the name of the combined dark
        masterdark.write(path / dark_file_name, overwrite = True )  # save the combined dark in path 
        
        if distribution == True: # compare the distribution of 1 dark and the combined one
        
            print ('PLOTTING BIN DATA COULD TAKE SOME TIME.')
        
            im1 = CCDData.read(calibrated_darks[0],unit='adu').data.flatten()
            im2 = masterdark.data.flatten()   # we read the data using CCDData and flatten() to be 1-Dimensional
        
            m = np.mean(im1)  # for setting the limit of the plot
            s = np.std(im1)
        
            lower = m - 0.1*s
            upper = m + 0.1*s
        
            plt.figure(figsize=(9, 9))
            hist(im1, bins=80, label='One single dark, ' + str(dark_times[0])+' [s]', alpha=1, histtype='step', range = (-100, 100)) # histogram of the data
            hist(im2, bins=80, label='{} Dark frames combined'.format(len(calibrated_darks))+ ' of ' + str(dark_times[0])+ ' [s]', alpha=1, histtype='step', range = (-100,100))
            plt.legend(loc=2, fontsize = 10)
            plt.grid()
            #plt.xlim(lower,upper)	
            plt.xlabel('Count level in image',fontsize = 20)
            plt.ylabel('Number of pixels with that count', fontsize = 20)
            plt.title(title, fontsize = 20)
            plt.savefig( path1 + 'Hist_Dark_Counts.png', dpi = 600,  bbox_inches='tight');
    
    if same_time == False:
        
        calibrated_darks = im_collection.files_filtered(object='dark', include_path=True, combined =None)

        masterdark = ccdp.combine(calibrated_darks,
                             method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6
                            )

        masterdark.meta['combined'] = True
    
        dark_file_name = 'masterdark.fit'

        masdarkter.write(path / dark_file_name, overwrite = True)
        
        if distribution == True: # compare the distribution of 1 dark and the combined one
        
            print ('PLOTTING BIN DATA COULD TAKE SOME TIME.')
        
            im1 = CCDData.read(calibrated_darks[0]).data.flatten()
            im2 = masterdark.data.flatten()   # we read the data using CCDData and flatten() to be 1-Dimensional
        
            m = np.mean(im1)  # for setting the limit of the plot
            s = np.std(im1)
        
            lower = m - 4*s
            upper = m + 4.5*s
        
            plt.figure(figsize=(9, 9))
            hist(im1, bins=5000, label='One single dark', alpha=0.5, range = (0,25)) # histogram of the data
            hist(im2, bins=5000, label='{} Dark frames combined'.format(len(calibrated_darks)), alpha=0.5, range = (0,25))
            plt.legend(fontsize = 13)
            plt.grid()
            plt.xlim(lower,upper)
            plt.xlabel('Count level in image',fontsize = 20)
            plt.ylabel('Number of pixels with that count', fontsize = 20)
            plt.savefig( path1 + 'Hist_Dark_Counts.png', dpi = 600,  bbox_inches='tight');
        
    
    #########################################################################################
    
    if compare == True:  # if true, will plot a single dark and the combined dark to compare them
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        show_image(CCDData.read(calibrated_darks[0],unit='adu').data, cmap='gray', ax=ax1, fig=fig)
        ax1.set_title('Single calibrated dark', fontsize = 20)
        show_image(masterdark.data, cmap='gray', ax=ax2, fig=fig)
        ax2.set_title('{} dark images combined'.format(len(calibrated_darks)), fontsize = 20)
       
       
       





def get_calibrated_flat(path_flat, path_master_dark, path_master_bias = None , same_time = True , folder = False, compare = False, cmap = 'gray' , warning = 'ignore' ):
    
    # This function will calibrated our flat frames by subtracting bias and dark frames.
    # we cannot just combined the flat frames, we need to calibrate them first, just as the scientific images.
    
    '''
       path_flat: path where your uncalibrated flat frames are. Make sure they have the same exposure time.
       
       folder: if you set it True it will show what is inside your folder, by default it is False.
       
       warning: decide you want to ignore warning or not. By default it is ignore. If you want to
                 see them, set warning = 'default'.
       
       same_time: If True, It will only calibrated flat frames by subtracting the master_dark of the same exposure time, here you dont need to subtract bis frame.
                  If False, It will calibrate the flats frames by subtractiong both master_dark and master_bias regardless the expossure time.
                  By default True.
       
       path_master_dark: path to the master dark obtained by combining (calibrated if needed) dark frames.
       
       path_master_bias: path to the master bias obtained by combining bias frames.
       
       compare: if True, you will compare the mask with the flat you use to get it. By default it is False.
        
       cmap: you can choose the colormap of the plot.
                  
       '''
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    #######################################################################################
    
    # here we modify the directory where the flat frames are so we create a new folder for the calibrated flats
    
    import shutil
    
    folder_path = Path(path_flat + 'calibrated_flats/')

    # Delete the existing folder and all its contents
    if folder_path.is_dir():
    	shutil.rmtree(folder_path)
    
    
    
    if path_flat[-1] == '/':                
        
        path_flat_new = path_flat + 'calibrated_flats/'
        path_calibrated_flats = Path(path_flat_new)
        path_calibrated_flats.mkdir()
        print('THE CALIBRATED FLATS ARE PLACE IN:   ' + path_flat_new)
    
    if path_flat[-1] != '/':
        
        path_flat_new = path_flat + '/calibrated_flats/'
        path_calibrated_flats = Path(path_flat_new)
        path_calibrated_flats.mkdir()
        print('THE CALIBRATED FLATS ARE PLACE IN:   ' + path_flat_new)
       
    
    #######################################################################################    
   
    path_flat = Path(path_flat)
    
    flats = ImageFileCollection(path_flat)
    
    path_master_dark = Path(path_master_dark)
    
    master_dark = CCDData.read(path_master_dark)
    
    if path_master_bias != None:
        
        path_master_bias = Path(path_master_bias)
    
        master_bias = CCDData.read(path_master_bias)
   
    #######################################################################################

    if folder == True:
        
        # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(flats.summary)   # display the files in the notebook
        
    ########################################################################################
    
    if same_time == True:
        
    
        flat_times = set(flats.summary['exptime']) # the exposure times of your dark frames
    
        flat_times_list = list(flat_times)
    
        if len(flat_times_list) > 1:  # we are checking if there are different exposure times which should not be the case
        
            raise Exception('YOU HAVE FLAT WITH DIFFERENT EXPOSURE TIME, PLEASE CHECK IT.')
            print(f'EXPOSURE TIME OF YOUR FLAT FRAMES ARE {flat_times_list} SECONDS.')
            #exit()  
        
        if len(flat_times_list) == 1:
    
            print(f'EXPOSURE TIME OF YOUR FLAT FRAMES IS {flat_times_list[0]} SECONDS.')
        
    #########################################################################################
    
    
        exp_time_master_dark = master_dark.header['exptime']
    
        if exp_time_master_dark != flat_times_list[0]: # we are checking if flat and dark have the same exposure time.
        
            raise Exception('MASTER DARK AND FLATS FRAMES HAVE DIFFERENT EXPOSURE TIME, PLEASE CHECK IT. THEY SHOULD HAVE THE SAME EXPOSURE TIME.')
            #exit()
        
    #########################################################################################
    
        for ccd, file_name in flats.ccds( object = 'flat', exptime=flat_times_list[0] ,     
                                          ccd_kwargs={'unit': 'adu'} ,   # CCDData requires a unit for the image if 
                                          return_fname = True ):         # Provide the file name too.                         
                                                                                        
            # Subtract the dark current 
        
            ccd = ccdp.subtract_bias(ccd, master_bias)
            
            ccd = ccdp.subtract_dark(ccd, master_dark,
                                 exposure_time='exptime', exposure_unit=u.second)

            # Save the result; there are some duplicate file names so pre-pend "calibrated_"
        
            ccd.write( path_calibrated_flats / ('calibrated_' + file_name))
        
    ##########################################################################################
    
    # here we plot the calibrated and ucalibrated flat. I have chosen it to be the first one because none of them should be more special than the other ones.
    
        if compare == True:
            
                flat_calib = ImageFileCollection(path_calibrated_flats)
            
                im2 = flat_calib.files_filtered(exptime = sorted(flat_times)[0], include_path = True)
            
                im1 = flats.files_filtered(exptime = sorted(flat_times)[0] ,include_path=True)
            
                no_calibrated = CCDData.read(im1[0],unit='adu')
            
                calibrated = CCDData.read(im2[0],unit='adu')
            
                fig, axes = plt.subplots(1, 2, figsize=(20, 10))

                show_image(no_calibrated, cmap=cmap, fig=fig, ax=axes[0])
                axes[0].set_title('Uncalibrated Flat', fontsize = 20)

                show_image(calibrated, cmap= cmap, fig=fig, ax=axes[1])
                axes[1].set_title('Calibrated Flat', fontsize = 20)
    
    #################################################################################################
                
    if same_time == False:
        
        for ccd, file_name in flats.ccds(object ='flat', return_fname=True, ccd_kwargs={'unit': 'adu'}):
        
            # Subtract the bias
            ccd = ccdp.subtract_bias(ccd, master_bias)
    
            # Subtract the dark current 
            ccd = ccdp.subtract_dark(ccd, master_dark, exposure_time='exptime', exposure_unit=u.second)

            # Save the result
            ccd.write( path_calibrated_flats / ('calibrated_' + file_name))
     
     ###################################################################################################
        
        if compare == True:
            
                flat_calib = ImageFileCollection(path_calibrated_flats)
            
                im2 = flat_calib.files_filtered(imagetyp = 'LIGHT', include_path = True)
            
                im1 = flats.files_filtered(imagetyp = 'LIGHT' ,include_path=True)
            
                no_calibrated = CCDData.read(im1[0],unit='adu')
            
                calibrated = CCDData.read(im2[0],unit='adu')
            
                fig, axes = plt.subplots(1, 2, figsize=(20, 10))

                show_image(no_calibrated, cmap=cmap, fig=fig, ax=axes[0])
                axes[0].set_title('Uncalibrated Flat')

                show_image(calibrated, cmap= cmap, fig=fig, ax=axes[1])
                axes[1].set_title('Calibrated Flat')











def get_master_flat(path, folder = False, same_time = True ,compare = False, distribution = False, cmap = 'gray' ,warning = 'ignore'):
    
    
    # This function generate a master dark frame by combining dark frames of the same exposure time.
    
    ''' path: path to the folder where the calibrated flat frames are located.
        Example:  path = '/home/felipe/Hiwi/Pycode/folder_flat'
        THE CODE WONT WORK WITH FILE.TIF FORMAT, WORK WITH FILE.FIT FORMAT.                                '''
    
    ''' folder = if you set it True it will show what is inside your folder, by default it is False.       '''
    
    ''' same_time = If True It will combined calibrated flat frames of the same exposure time.
                    If False it will combined all the calibrated flat frames regardless the exposure time.
                    By default True.'''
    
    ''' compare = if you set it True it will plot a single calibrated flat and the combined flat to compare them.
                  By default it is False.       
                                                                             '''
    ''' warning = decide you want to ignore warning or not. By default it is ignore. If you want to
                  see them, set warning = 'default'  
                                                                        '''
    ''' so far, for this code, the combine files cannot be already in the folder, otherwise errors will
        arise. For that, if there is a file with the HEADER = 'combined', the code will stop.  
        Make sure, a combined dark file is not in the folder.                                              '''
   
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    path1 = path
    path = Path(path) 
    
    cflats = ImageFileCollection(path)   # read all the files in the folder given by path
    
    ##########################################################################################
    
    #for im in cflats.hdus(combined=True):  # seek if the combined flat is already there,
     #   if im != None:
      #      print('COMBINED FLAT IS ALREADY IN THE FOLDER!')
       #     exit()                                 # if it is already there, quit the code
         
    #########################################################################################
    
    if folder == True:          # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(cflats.summary)   # display the files in the notebook
    
    flat_times = set(cflats.summary['exptime']) # the exposure times of your flat frames
    
    ##########################################################################################
    
    def inv_median(a):   # we defined this function to normalize our flats
        return 1 / np.median(a)
    
    #########################################################################################
    
    if same_time == True:
        
        for exp_time in sorted(flat_times):
        
            calibrated_flats = cflats.files_filtered(object = 'flat', exptime=exp_time, include_path=True, combined = None)
            combined_flat = ccdp.combine(calibrated_flats, method='median', mem_limit=350e6)


            combined_flat.meta['combined'] = True
            flat_file_name = 'combined_flat_{:6.3f}.fit'.format(exp_time)  # set the name of the combined dark
            combined_flat.write(path / flat_file_name, overwrite = True)
            warnings.warn("We will overwrite a new combined flat if it was already there!")
            
        if distribution == True: # compare how the median and mean change in each image, to make clear why we need to normalize
        
            print('PLOTTING BIN DATA COULD TAKE SOME TIME.')
         
            median_count = [np.median(data) for data in cflats.data(exptime = sorted(flat_times)[0])]
            mean_count = [np.mean(data) for data in cflats.data(exptime = sorted(flat_times)[0])]
            med = list(map('{:.0f}'.format,mean_count))
            plt.plot(median_count, label='median')
            plt.plot(mean_count, label='mean')
            plt.xlabel('Image number', fontsize  =20)
            plt.ylabel('Count (ADU)', fontsize  =20)
            plt.title('Pixel value in calibrated flat frames', fontsize  =20)
            plt.legend(loc = 'best', fontsize  =14)
            print(f'THE MEDIAN OF EACH IMAGE ARE {med}')
            plt.xlim(0, len(median_count)-1)
            plt.savefig( path1 + 'pixel_value_calibrated_flats.png', dpi = 600,  bbox_inches='tight');
    
    if same_time == False:
    
        calibrated_flats = cflats.files_filtered(object = 'flat', combined = None, include_path=True)
        combined_flat = ccdp.combine(calibrated_flats,
                                     method='average', scale=inv_median,
                                     sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,
                                     mem_limit=350e6
                                    )

        combined_flat.meta['combined'] = True
        flat_file_name = 'combined_flat.fit'  # set the name of the combined dark
        combined_flat.write(path / flat_file_name, overwrite = True)
        warnings.warn("We will overwrite a new combined flat if it was already there!")
        
        if distribution == True: # compare how the median and mean change in each image, to make clear why we need to normalize
        
            print ('PLOTTING BIN DATA COULD TAKE SOME TIME.')
         
            median_count = [np.median(data) for data in cflats.data(combined = None)]
            mean_count = [np.mean(data) for data in cflats.data(combined = None)]
            med = list(map('{:.0f}'.format,mean_count))
            plt.plot(median_count, label='median')
            plt.plot(mean_count, label='mean')
            plt.xlabel('Image number', fontsize  =15)
            plt.ylabel('Count (ADU)', fontsize  =15)
            plt.title('Pixel value in calibrated flat frames', fontsize  =15)
            plt.legend(fontsize  =15)
            print(f'THE MEDIAN OF EACH IMAGE ARE {med}')
            plt.savefig( path1 + 'pixel_value_calibrated_flats.png', dpi = 600,  bbox_inches='tight');
        
        
    #########################################################################################
    
    if compare == True:  # if true, will plot a single calibrated flat and the combined normalized flat to compare them
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        show_image(CCDData.read(calibrated_flats[0],unit='adu').data, cmap=cmap, ax=ax1, fig=fig)
        ax1.set_title('Single calibrated flat', fontsize = 20)
        show_image(combined_flat.data*inv_median(combined_flat.data), cmap=cmap, ax=ax2, fig=fig)
        ax2.set_title('{} flat images combined and normalized'.format(len(calibrated_flats)), fontsize = 20)
    










def get_hot_pixels(path_im, path_fol, sigma = 5, cenfunc = 'median', stdfunc = 'mad_std' ,compare = False, hot_pixels = False ,cmap = 'gray' ,warning = 'ignore'):
    
    # This function generate a mask for hot pixels. Use a dark frame of long exposure time.
    
    ''' path_im: path to where the dark frame is located.
        Example:  path = '/home/felipe/Hiwi/Pycode/folder_dark'
        THE CODE WONT WORK WITH FILE.TIF FORMAT, WORK WITH FILE.FIT FORMAT.    
                                     '''
    ''' path_fold: path to where the mask should be saved. Highly recomended to save in the save folder as
                   the dark frames.    
                                                                                        '''
    ''' sigma: the number of standard deviation you would like to clip the data. By default 5 sigmas.      
     '''
    ''' cenfunc: The statistic or callable function/object used to compute the center value for the clipping.
                 you can choose {'median’, ‘mean'}.   
                                                                       '''
    ''' stdfunc: The statistic or callable function/object used to compute the standard deviation about 
                the center value. You can choose {‘std’, ‘mad_std’}.       
                                                 '''
    ''' cmap: to change colormap colors. By default 'gray'. A good option is 'viridis'. 
                        '''
    ''' warning = decide you want to ignore warning or not. By default it is ignore. If you want to
                  see them, set warning = 'default'                                                        '''
    
    
    #############################################################################################
    
    warnings.simplefilter( warning , category=AstropyWarning) # for avoiding the warnings from astropy
    
    path_im = Path(path_im)
    
    path_fol = Path(path_fol)
    
    #############################################################################################
    
    #im_collection = ImageFileCollection(path_fol)  # creating a collection of images from all the files
                                                   # in path_fol
    #for im in im_collection.hdus(imagetyp= 'dark mask'): # if the mask is already there, quit the code
    #    if im != None:
    #        print('DARK MASK IS ALREADY IN THE FOLDER!')
    #        exit()
    
    #############################################################################################
    
    im = CCDData.read(path_im) # read the image.fit
    
    im_clipped = sigma_clip(im , sigma = sigma , cenfunc = cenfunc, stdfunc=stdfunc) # we use sigma_clip
    # from astropy to calculate the hot pixels. Those pixels with a given amount of stardard deviations
    # away from the mean are marked as hot pixels and then we create a mask using them.
    
    #show_image(im_clipped.mask,  cmap=cmap , is_mask=True) # we show the mask
    
    ##############################################################################################
    
    data1 = im_clipped.mask # this is for saving the mask we have just created
    
    if hot_pixels == True:
        
        print(f'WE HAVE FOUN {data1.flatten().sum()} HOT PIXELS.')
        
    
    mask_as_ccd = CCDData(data1.astype('uint8'), unit=u.dimensionless_unscaled) # saving the mask
    
    mask_as_ccd.header['imagetyp'] = 'dark mask'  # adding the Header 'imagetyp' = 'dark mask'
    
    mask_as_ccd.write( path_fol / 'mask_from_dark_current.fits' , overwrite=True)  # saving in path_fol
    
    ###############################################################################################
    
    if compare == True:
        
        fig, axes = plt.subplots(1, 2, figsize=(15, 10))

        show_image(im, cmap= cmap, fig=fig, ax=axes[0])
        axes[0].set_title('Single Dark Frame', fontsize = 20)

        show_image(data1, cmap = cmap, fig=fig, ax=axes[1], is_mask=True)
        axes[1].set_title('Hot Pixels Mask', fontsize = 20);
    
    
    



        
def get_bad_pixels(path1, path_fol , path2 = None , compare = False, cmap = 'gray', bad_pixels = False, warning = 'ignore'):
    
    # This function will calculate the BAD PIXELS in the CCD camera.
    
    ''' This function will create a mask for bad pixels. It uses flat frames to calculate it.
        There are two methods, the first one is using only single flat and the secod one is using two flats.
        If you have flats with same average counts, then it is better to use only one flat, if you have
        two flats with different counts, you can use the second method. By default if you only give one flat
        path, it will use the method using only one flat.
        
        path1: path to the first flat. Mandatory.
        
        path_fol: path to the folder where you want to save it. Mandatory.
        
        path2: path to the second flat. Optional.
        
        compare: if True, you will compare the mask with the flat you use to get it. By default it is False.
        
        cmap: you can choose the colormap of the plot.
        
        bad_pixels: if it is True, it will print the amount of bad pixels found. By default it is False.
        '''
    
    ##########################################################################################################
    
    warnings.simplefilter( warning , category=AstropyWarning) # for avoiding the warnings from astropy
    
    path1 = Path(path1)
    
    path_fol = Path(path_fol)
    
    ##########################################################################################################
    
    if path2 == None:
        
        print('IT COULD TAKE SOME TIME, BUT YOU JUST NEED TO DO IT ONCE.')
        
        flat = CCDData.read(path1)  # we are reading the chosen flat frame.
        mask_bad_pixels = ccdp.ccdmask(flat)  # we calculate the bad pixels by using the function ccdmask from ccdp python module.
        
        mask_as_ccd = CCDData(data=mask_bad_pixels.astype('uint8'), unit=u.dimensionless_unscaled)  # we are saving our mask to be use later.
        mask_as_ccd.header['imagetyp'] = 'flat mask'   # we add the header imagetyp = flat mask to specify that this is the mask from flat frames.
        mask_as_ccd.write( path_fol / 'mask_from_ccdmask.fits' , overwrite=True)
        
        if bad_pixels == True:
            
            print(f'WE HAVE FOUN {mask_bad_pixels.flatten().sum()} BAD PIXELS.')
        
        if compare == True:
            
            fig, axes = plt.subplots(1, 2, figsize=(15, 10))

            show_image(flat, cmap= cmap, fig=fig, ax=axes[0])
            axes[0].set_title('Single calibrated flat')

            show_image(mask_bad_pixels, cmap = cmap, fig=fig, ax=axes[1], is_mask=False)
            axes[1].set_title('Derived mask');
    
    ##########################################################################################################
    
    if path2 != None:  # we start by using the method which compare two flat flats, we check if a path for a second flat frame was given.
        
        print('IT COULD TAKE SOME TIME, BUT YOU JUST NEED TO DO IT ONCE.')
        
        path2 = Path(path2)
        flat1 = CCDData.read(path1)
        flat2 = CCDData.read(path2)
        
        ratio = flat2.divide(flat1)  # the method used by cccdmask function requires that we need to divide the two flat frames.
        
        #print(ratio)
        
        print(f'THE RATIO IS {ratio.data.mean()}')
        
        mask_bad_pixels = ccdp.ccdmask(ratio)  # we calculate the bad pixels using ccdmask from the ratio.
        
        mask_as_ccd = CCDData(data=mask_bad_pixels.astype('uint8'), unit=u.dimensionless_unscaled)   # saving the mask.
        mask_as_ccd.header['imagetyp'] = 'flat mask'
        mask_as_ccd.write( path_fol / 'mask_from_ccdmask.fits' , overwrite=True)
        
        if bad_pixels == True:  # calculating and printing the number of bad pixels.
            
            print(f'WE HAVE FOUN {mask_bad_pixels.flatten().sum()} BAD PIXELS.')
        
        if compare == True:
            
            fig, axes = plt.subplots(1, 2, figsize=(20, 10))   # comparing the resulting mask and the ratio image

            show_image(ratio, cmap=cmap, fig=fig, ax=axes[0], show_colorbar=False)
            axes[0].set_title('Ratio of two flats', fontsize = 20)

            show_image(mask_as_ccd, cmap= cmap, fig=fig, ax=axes[1], show_colorbar=False, percl=99.95, is_mask = True)
            axes[1].set_title('Bad Pixels Mask', fontsize = 20)
        
        
        
        
        
        
                

def get_combined_mask(darkmask, flatmask, path_fol , useless_pixels = False , show_mask = False, cmap = 'viridis'):
    
    # This function will combined the mask obtained by finding the hot pixels and bad pixels.
    '''
       darkmask: path to the mask obtained from exceed of dark current, hot pixels. Mandatory.
       
       flatmask: path to the mask obtained from flat frames, bad pixels. Mandatory.
       
       path_fol: path to the directory where you want to save the mask. Mandatory.
       
       useless_pixels: if True, will print the total amount of pixel that will be neglected. By default it is False.
       
       show_mask: if True, will plot the new mask. By default it is False.
       
       cmap: you can change colormap of the plot.
       
       '''
    
    darkmask = Path(darkmask)
    flatmask = Path(flatmask)
    path_fol = Path(path_fol)
    
    ############################################################################################################
    
    mask_ccdmask = CCDData.read(flatmask, unit=u.dimensionless_unscaled)
    mask_ccdmask.data = mask_ccdmask.data.astype('bool')

    mask_hot_pix = CCDData.read(darkmask, unit=u.dimensionless_unscaled)
    mask_hot_pix.data = mask_hot_pix.data.astype('bool')
    
    #############################################################################################################
    
    combined_mask = mask_ccdmask.data | mask_hot_pix.data
    
    useless_pixel = combined_mask.sum()
    
    total_pixels = combined_mask.shape[0]*combined_mask.shape[1]
    
    percentage_useless_pixels = (useless_pixel/total_pixels)*100
    
    #############################################################################################################
    
    if useless_pixels == True:
        
        print(f'WE HAVE FOUND IN TOTAL {useless_pixel} USELESS PIXELS, WE ARE MASKING ROUGHLY {percentage_useless_pixels:.2f}% OF THE PIXELS.')
    
    if show_mask == True:
        
        show_image(combined_mask,  cmap=cmap , is_mask = True) # we show the mask
        plt.title(f'Combined Mask, {percentage_useless_pixels:.2f}% Pixels Not Working ',fontsize = 19)
        
    #############################################################################################################
    
    combined_mask = CCDData(data=combined_mask.astype('uint8'), unit=u.dimensionless_unscaled)
    combined_mask.header['imagetyp'] = 'combined mask'
    combined_mask.write( path_fol / 'combined_mask.fits' , overwrite=True)
    
    
    
    




def flatness(path, size, position, num_bins, x1_min, x1_max, x2_min, x2_max, cmap = 'gray', boxcolor = 'red'):

    # This function will measure how flat are the flatframes. It only evaluate one flat per time but can be extended to do it over a complete folder of flatframes.
    # In principle I would be ideal to do it with the master flat frames. 
    
    flat = CCDData.read(Path(path), unit = 'adu').data  # reading a single flat image
    flat_data = flat.flatten()  # flattening to have a 1D array instead of 2D arrays. Make easier to do histograms and for analysing the data.
    
    si = u.Quantity(size, u.pixel)  # defining the size of the square 
        
    cutout_zone = Cutout2D(flat, position, si)   # we use the function Cutout2D from ccdproc to cut the image
    
    cutout_data = cutout_zone.data.flatten()  # getting the data, Cutout2D return a more complex object than just an array
     
    #############################################################################################################################
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))     # from here we plot the flat together with the cutout region to visualize the cutout zone.
    cutout_zone.plot_on_original(ax = ax1, color= boxcolor)

    show_image(flat, cmap= cmap, ax=ax1, fig=fig)
    ax1.set_title('Flat Frame', fontsize = 20)
    
    show_image(cutout_zone.data, cmap= cmap, ax=ax2, fig=fig)
    ax2.set_title('Cutout Region from Flat Flame', fontsize = 20)  
    
    ##############################################################################################################################
    
    
    def gaussian(x, amplitude, mean, stddev):
        return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)   # definition of a gaussian distribution function for fittig
    
    ##############################################################################################################################
    
    bin_edges1 = np.linspace(x1_min, x1_max, num_bins + 1)      # defining the bin edges of the histograms
    bin_centers1 = (bin_edges1[1:] + bin_edges1[:-1]) / 2       # defining the bin centers for the histograms
    counts1, _ = np.histogram(flat_data, bins=bin_edges1, density=True)  # counts per bin
    indices1 = np.where((bin_centers1 >= x1_min) & (bin_centers1 <= x1_max))[0]  # indecing
    mean1 = np.sum(bin_centers1[indices1] * counts1[indices1]) / np.sum(counts1[indices1])  # calculting the mean of the selected region
    stddev1 = np.sqrt(np.sum((bin_centers1[indices1] - mean1) ** 2 * counts1[indices1]) / np.sum(counts1[indices1]))  # calculating standard deviation of the selected region
    initial_params1 = [np.max(counts1), mean1, stddev1]   # we use those values for mean and sigma as expected values for the fitting function
    popt1, pcov1 = curve_fit(gaussian, bin_centers1[indices1], counts1[indices1], p0=initial_params1)  # getting best fitting values from curve_fit function 
    
    #############################################################################################################################
    
    bin_edges2 = np.linspace(x2_min, x2_max, num_bins + 1)
    bin_centers2 = (bin_edges2[1:] + bin_edges2[:-1]) / 2
    counts2, _ = np.histogram(cutout_data, bins=bin_edges2, density=True)
    indices2 = np.where((bin_centers2 >= x2_min) & (bin_centers2 <= x2_max))[0]
    mean2 = np.sum(bin_centers2[indices2] * counts2[indices2]) / np.sum(counts2[indices2])
    stddev2 = np.sqrt(np.sum((bin_centers2[indices2] - mean2) ** 2 * counts2[indices2]) / np.sum(counts2[indices2]))
    initial_params2 = [np.max(counts2), mean2, stddev2]
    popt2, pcov2 = curve_fit(gaussian, bin_centers2[indices2], counts2[indices2], p0=initial_params2)
    
    ##############################################################################################################################
    
    fit_info1 = [f'$\\mu = {popt1[1]:.0f}$', f'$\\sigma = {popt1[2]:.0f}$']    # defining best fitting values to place them in the plots
    fit_info2 = [f'$\\mu$ = {popt2[1]:.0f}', f'$\\sigma$ = {popt2[2]:.0f}']
    
    #############################################################################################################################
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))    # plotting the histograms for the whole data and just for the selected region to compare how the fitting is improved
      
    ax2.hist(cutout_data, bins=bin_edges2, alpha=0.6, density=True, color = "skyblue", edgecolor='black')
    ax2.plot(bin_centers2[indices2], gaussian(bin_centers2[indices2], *popt2), 'r-', label='fit')
    ax2.legend(loc='best', title="\n".join(fit_info2), fontsize  = 15, title_fontsize=15)
    ax2.set_title('Count Distribution Cutout Zone', fontsize = 20)
    ax2.set_xlabel('Count Level', fontsize = 20)
    ax2.set_yticks([])
    
    ###############################################################################################################################
    
    ax1.hist(flat_data, bins=bin_edges1, alpha=0.6, density=True, color = "skyblue", edgecolor='black')
    ax1.plot(bin_centers1[indices1], gaussian(bin_centers1[indices1], *popt1), 'r-', label='fit')
    ax1.legend(loc='best', title="\n".join(fit_info1), fontsize  = 15, title_fontsize=15)
    ax1.set_title('Count Distribution Full Flat Frame', fontsize = 20)
    ax1.set_xlabel('Count Level', fontsize = 20)
    ax1.set_yticks([])


    # Calculate the errors on the fitted parameters
    perr2 = np.sqrt(np.diag(pcov2))
    
    #print(f'$\\mu = {popt1[1]:.0f}  \\pm {popt1[2]:.0f}$')
    
    #return tuple(param_strings2)
    






def get_calibrated_ccd_image(path_science, path_master_dark, path_master_flat, matrix_error, path_mask, path_master_bias ,gain = None , compare = False,  readnoise = None, cmap = 'hot', warning = 'ignore'):
    
    # we apply all the calibration images to the measurements so we obtained a reduced image which we want extract phyics from.
    
    '''
    
    path_science: path to the image(s) you need to calibrate.
    
    path_master_bias: path to the master bias, only needed if you are working with images of different exposure time.
                      Therefore you can subtract the bias level using dark frames. By default None.
    
    path_master_dark: path to the master dark, which you should obtain by using the function get_master_dark.
    
    path_master_flat: path to the master flat, which you should obtain by using the function get_master_flat.
    
    path_mask: path to the mask which contains all the pixels that should not be taken into account for analysis. You shoud get it by using get_combined_mask function.
    
    gain: a gain value will multiply the image so that you scale the picture. By default is None.
    
    readnoise: you can specify the readnoise of your camera, it could maybe change in some conditions. By default None.
    
    compare: if True, it will plot the image before and after the reduction, so you can compare them. By default it is False. 
    
    cmap: you can change the colormap. By default it is viridis.
    
    distribution: if True, will show the distribution (a histogram) of the image after and before the calibration. By default it is False.
    
    '''
    
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    ###########################################################################################
    
    # here we modify the directory where the scientific pictures are so we create a new folder for the calibrated pictures
    
    import shutil
    
    folder_path = Path(path_science + 'calibrated_science/')

    # Delete the existing folder and all its contents
    if folder_path.is_dir():
        shutil.rmtree(folder_path)
    
    
    if path_science[-1] == '/':                
        
        path_science_new = path_science + 'calibrated_science/'
        path_calibrated_science = Path(path_science_new)
        path_calibrated_science.mkdir()
        print('THE CALIBRATED SCIENTIFIC PICTURES ARE PLACE IN:   ' + path_science_new)
    
    if path_science[-1] != '/':
        
        path_science_new = path_science + '/calibrated_science/'
        path_calibrated_science = Path(path_science_new)
        path_calibrated_science.mkdir()
        print('THE CALIBRATED SCIENTIFIC PICTURES ARE PLACE IN:   ' + path_science_new)
    
    ###########################################################################################
    
      
    path_science = Path(path_science) 
    path_master_dark = Path(path_master_dark)
    path_master_flat = Path(path_master_flat)
    path_master_bias = Path(path_master_bias)
    path_mask = Path(path_mask)
    
    images = ImageFileCollection(path_science)     # reading the data we already obtained from previus processes.
    master_dark = CCDData.read(path_master_dark)
    master_flat = CCDData.read(path_master_flat)
    mask = CCDData.read(path_mask, unit = 'adu')
    mask_array = mask.data
    master_bias = CCDData.read(path_master_bias)
        
    def uncertainty(List):   
        l = [np.square(m) for m in List]
        s = sum(l)
        return np.sqrt(s)
    
    matrix_error = uncertainty(matrix_error)
        
    for science, file_name in images.ccds(object = 'science', reduced = None, combined = None ,return_fname=True, ccd_kwargs=dict(unit='adu')):
    
        science.mask = mask_array
        
        reduced_image = ccdp.subtract_bias(science, master_bias)
        
        reduced_image = ccdp.subtract_dark(reduced_image, master_dark, exposure_time='exptime', exposure_unit=u.second)
            
        reduced_image = ccdp.flat_correct(reduced_image, master_flat, min_value = 0.9)
        
        reduced_image.uncertainty = matrix_error
        
        if gain != None:
            
            reduced_image = ccdp.gain_correct(reduced_image , gain*u.electron/u.adu)
            
        reduced_image.meta['reduced'] = True  # add to the reduced file the HEADER combined = 'True'
    
        reduced_image.write( path_calibrated_science / ('reduced_' + file_name), overwrite = True )  # save the reduced image in path_cdd_image 
        
        if compare == True:  
        
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
            show_image(science.data, cmap=cmap, ax=ax1, fig=fig)
            ax1.set_title('Uncalibrated Image', fontsize = 20)
            show_image(reduced_image.data, cmap=cmap, ax=ax2, fig=fig)
            ax2.set_title('Calibrated Image', fontsize = 20)








def get_flats_for_test(path, compare = False, folder = False ,cmap = 'gray' , warning = 'ignore'):
    
    # this function will calibrate and prepare the flats frames for calculating the gain using the function get_gain_factor
    # since you will need to take different picture with different exposure time, it will be more convinient to write a function
    # that will take the folder with all these pictures and processes them at once, instead of using get_calibrated_flat function for each pair
    
    '''
       path: path to the folder that should contains the pair of flats and the master dark frames with the same exposure time.
       
       folder: if you set it True it will show what is inside your folder, by default it is False.
       
       warning: decide you want to ignore warning or not. By default it is ignore. If you want to
                 see them, set warning = 'default'
       
       compare: if True, you will display the flats before and after calibration. By default it is False.
       
       cmap: you can choose the colormap of the plot.
    '''
    
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    ########################################################################################
    
    import shutil
    
    folder_path = Path(path + 'calibrated_flats/')

    # Delete the existing folder and all its contents
    if folder_path.is_dir():
    	shutil.rmtree(folder_path)
    
    
    
    if path[-1] == '/':                 # for creating a new folder, we need to create the a new directory
         
        path_flat_new = path + 'calibrated_flats/'
        path_calibrated_flats = Path(path_flat_new)
        path_calibrated_flats.mkdir()
        print('THE CALIBRATED FLATS ARE PLACE IN:   ' + path_flat_new)
    
    if path[-1] != '/':                # we are checking the structure of the path so we can create the new directory
        
        path_flat_new = path + '/calibrated_flats/'
        path_calibrated_flats = Path(path_flat_new)
        path_calibrated_flats.mkdir()
        print('THE CALIBRATED FLATS ARE PLACE IN:   ' + path_flat_new)
    
    #########################################################################################
    
    path = Path(path)
    
    im_collection  = ImageFileCollection(path)   # we create an image collection to work with all the files at once
   
    #########################################################################################

    for im in im_collection.hdus(SUBTRACT_DARK = 'subdark'):  # seek if the calibrated flats are already there,
        if im != None:
            raise Exception('CALIBRATED FLAT IS ALREADY IN THE FOLDER')
    
    ##########################################################################################
    
    if folder == True:
        
        # to show the files in the folder if you set folder = True
        print('DATA IN THE FOLDER')
        display(im_collection.summary)   # display the files in the notebook
    
    ##########################################################################################
    
    flat_times = set(im_collection.summary['exptime']) # the exposure times of your flat frames
    
    print(f'THE EXPOSURE TIMES ARE {sorted(flat_times)} SECONDS.')
    
    ##########################################################################################
    
    combined_darks = {ccd.header['exptime']: ccd for ccd in im_collection.ccds(object = 'dark', combined=True)}  # make sure there are not other combined frames but dark frames
    
    masterbias = CCDData.read(Path(im_collection.files_filtered( object = 'bias', include_path=True, combined = True)[0]))
    #print(masterbias.data)
    
   
                                         
                                         
    
    for ccd, file_name in im_collection.ccds(object = 'flat', combined = None,         
                                         ccd_kwargs={'unit': 'adu'}, # CCDData requires a unit for the image if 
                                                                     # it is not in the header
                                         return_fname=True           # Provide the file name too.
                                        ):   
    # Subtract the dark current 
   
        
        
        ccd_calibrated = ccdp.subtract_bias(ccd, masterbias)
        
        ccd_calibrated = ccdp.subtract_dark(ccd, combined_darks[ccd.header['exptime']], exposure_time='exptime', exposure_unit=u.second)  # Subtract the dark and bias by using a master dark of the same 
                                                                                                                                          # exposure time
                                                                               # It makes sure that the time of the flats are the same of the master darks
         
        ccd_calibrated.write( path_calibrated_flats / ('calibrated_' + file_name) , overwrite = True)
        #print(ccd.header['exptime'])
    
        if compare == True:  
        
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
            show_image(ccd.data, cmap=cmap, ax=ax1, fig=fig)
            ax1.set_title('Uncalibrated Flat Frame', fontsize = 16)
            show_image(ccd_calibrated.data, cmap=cmap, ax=ax2, fig=fig)
            ax2.set_title('Calibrated Flat Frame + {}'.format(ccd.header['exptime']), fontsize = 16)
    







def get_gain_factor(path , position , size , gain , title = 'Variance diagram for CCD: slope yields gain' , compare = False , save = False ,cmap = 'gray' , warning = 'ignore'):
    
    # this function calculate the gain of your CCD camera by using flat frames taken for different exposure times.
    # We need to calculate the gain of the cameras since we need the cameras to convert the amount of photon to count by the same factor.
    # for using this, you have to have already a folder with all the flat frames with different exposure times calibrated.
    ''' 
       path: path to the folder where you have place the pair of flats. Please create and additional folder and place the flats there. The flats have to be already calibrated.
       
       position: coordinates where you think there is a zone with no bad or hot pixels. It have the following structure position = (x,y).
       
       size: size of the cutout. highly recommended to be an square with size between 40 to 100 pixels. Example: size = (y,x).
    
       gain: the gain given by the manufacturer. It will use for the chi2 calculation of the gain.
       
       cmap:  you can decide if you want to change the colormap. By default it is 'gray'.
       
       compare: will display the average flat frame with the zone you have cutout and compare it with the zone you have cutout. By default it is False.
       
       save: if True, it will save the Variance diagram for the ccd in path. By default it is False.
    
    '''
    
    warnings.simplefilter( warning , category=AstropyWarning)   # for avoiding the warnings
    
    #######################################################################################################################
    
    path1 = path
    
    path = Path(path)
    
    flats = ImageFileCollection(path)  # creating a collection with all the images inside the folder placed in path.
    
    flats_times = set(flats.summary['exptime']) # obtaining the exposure times. You need at least three.
    
    print(f'THE EXPOSURE TIMES ARE {sorted(flats_times)} SECONDS.')
    
    ########################################################################################################################
    
    if len(sorted(flats_times))<3:
        raise ValueError('YOU NEED AT LEAST 3 PAIR OF FLATS WITH DIFFERENT EXPOSURE TIME.')
    
    x = [] 
    y = []
    
    gain_m = gain
    ########################################################################################################################
    
    for exp_time in sorted(flats_times):
        
        flat_time = flats.files_filtered(object = 'flat', exptime=exp_time,include_path=True)  # We create a loop to run over images with the same exposure time
        
        if len(flat_time) > 2: # you need only two flat frames
            raise ValueError('THERE ARE MORE THAN 2 FLAT FRAMES WITH THE SAME EXPOSURE TIME.')
        
        if len(flat_time) < 2: # you need two flat frames
            raise ValueError('YOU NEED 2 FLAT FRAMES WITH THE SAME EXPOSURE TIME, YOU HAVE LESS OF 2.')
        
        flat1 = flat_time[0]  # reading the path from the list flat_time
        flat2 = flat_time[1]
        
        flat1 = CCDData.read(Path(flat1)).data  # we obtain the array data since CCDData is a more comple object
        flat2 = CCDData.read(Path(flat2)).data
        
        average = (flat1 + flat2)*0.5   # calculating the average and difference to be used in the method that will help us to find the gain
        
        difference = flat2 - flat1
        
        #show_image(average, cmap=cmap)
        
        #plt.title('Average Flat Frame', fontsize = 16)
        
        #plt.pause(0.001)
        
        #x = float(input('x COORDINATE FOR THE AREA WITH NO BAD PIXELS OR BAD PIXELS:  '))
        
        #y = float(input('y COORDINATE FOR THE AREA WITH NO BAD PIXELS OR BAD PIXELS:  '))
        
        po = position    # position of the good zone of the picture
        si = u.Quantity(size, u.pixel)  # pixels # size of the picture
        
        cutout_average = Cutout2D(average, po, si)   # we use the function Cutout2D from ccdproc to cut the image
        cutout_difference = Cutout2D(difference, po, si)
        
    #########################################################################################################################################
    
        if compare == True:  # in case you want to see what you have cut
        
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
            cutout_average.plot_on_original(ax = ax1, color='red')
            show_image(average, cmap=cmap, ax=ax1, fig=fig)
           
            
            ax1.set_title('Average Flat Frame', fontsize = 20)
            show_image(cutout_average.data, cmap=cmap, ax=ax2, fig=fig)
            ax2.set_title('Cutout Region for the Average Flat Flame', fontsize = 20)
        
        mean_level = np.mean(cutout_average.data)  # we caculate the mean_level to be used for calculating the gain
        variance = np.std(cutout_difference.data)*np.std(cutout_difference.data)/2
        
        x.append(mean_level)  # we put the data in a list 
        
        y.append(variance)
        
    ######################################################################################################################################
        
    x = np.array(x)  # converting the list into an array
    y = np.array(y)
    y_err = np.sqrt(y)/y  # we need error for our y data # poisson error
    
    def linear(x,a,b):  # model to be fitted. Here we use a linear function
        return a*x + b
    
    least_squares = LeastSquares(x, y, y_err,  model = linear) # we use LeastSquares function from iminuit numpy module for fitting
    
    m = Minuit(least_squares,a= 1/gain, b=0)  # starting values for a and b

    m.migrad()  # finds minimum of least_squares function
    m.hesse()   # accurately computes uncertainties

    
    # draw data and fitted line
    plt.figure(figsize = (10,8))
    plt.errorbar(x, y, y_err, fmt="o", label="data")
    plt.plot(x, linear(x, *m.values), label="fit")
    
    
    gain = m.values[0]**-1	

    diff = 100*(abs(gain_m - gain))/gain_m
    # display legend with some fit info
    fit_info = [f'Found gain $= {gain:.3f}$ [$e^-/adu$]', f'Manufacturer gain $= {gain_m:.2f}$ [$e^-$/adu]', f'We found a {diff:.2f}% difference']
    #for p, v, e in zip(m.parameters, m.values, m.errors):
    #    fit_info.append(f"{p} = ${v:.2f} \\pm {e:.2f}$")

    plt.title(title, fontsize=20)
    plt.ylabel('Variance inside box' , fontsize = 20)    
    plt.xlabel('Mean counts inside box', fontsize = 20)
    plt.legend(title="\n".join(fit_info), fontsize  = 15, title_fontsize=15)
    plt.xlim(min(x)-min(x)*0.02,max(x)+max(x)*0.01)
    
    ############################################################################################################################################
    
    # we calculate the gain from the plot, the gain is given by the inverse of the slope. gain = 1/a according to our model
    
    print(f'YOUR CAMERA HAS A GAIN OF {gain:.3f} electrons/ADU.')         
    
    ###########################################################################################################################################
   
    if save == True:
        
        plt.savefig( path1 + 'Variance_diagram_for_ccd.png', dpi = 600,  bbox_inches='tight');
        
        
        
        
        
        
        

def linearity_test_cdd(path, mask_path = None ,adu_max = None, adu_min = None, save = False ,warning = 'ignore'):
    
    # this code will provide information about the linearity of your ccd camera
    # you need to take flat frames with different exposure time, also the corresponding dark frames to calibrate the flat images
    # once they are calibrated you need to calculate the mean counts of the flats and then plot them vs time from which you should get a linear function
    # We do this by applying a chi2 test where the data error is considered as poisson error data_error = sqrt(data)/data
    # After the analysis is done, you will be able to see the range of counts in which you could work at
    
    '''
       Path:  path to the folder containing all the flat frames already calibrated. You only need one flat frame per time.
       
       adu_max: put the maximum number of counts you want to consider. If None, it will be approx 65000 counts.
       
       adu_min: put the minimum number of counts you want to consider. If None, it will zero the minimum amount of counts.
       
       save: if True, it will save the Linearity test plot in path. By default it is False.
    
       warning: decide you want to ignore warning or not. By default it is ignore. If you want to
                 see them, set warning = 'default'.
       
    '''
    
    warnings.simplefilter( warning , category=AstropyWarning)
    
    path1 = path 
    path = Path(path)
    
    calibrated_flats = ImageFileCollection(path)  # creating a collection of images from the data in the folder placed at path
    
    #######################################################################################################################
    
    mean_level = []   # empty list to save the mean_level of each flat frame
    
    data_err = []     # empty list to save the possion error of our data
    
    #######################################################################################################################
    
    flat_times = set(calibrated_flats.summary['exptime']) # the exposure times of your flat frames
    
    exptimes_flats = sorted(flat_times)
    
    print(f'YOUR EXPOSURE TIMES ARE {exptimes_flats} SECONDS.')
    
    #######################################################################################################################
    
    if mask_path == None:
    
        for exp_time in sorted(flat_times):
        
            # we filter them by the exposure time, this will ensure that we will use only a flat of a given time
        
            flat = calibrated_flats.files_filtered(object = 'flat', exptime=exp_time,include_path=True, combined=None)
        
            if len(flat) > 1:
                raise Exception('THERE ARE AT LEAST 2 FLAT FRAMES WITH THE SAME EXPOSURE TIME, YOU NEED ONLY ONE.')
            
            flat = CCDData.read(Path(flat[0])).data
        
            mean = np.mean(flat)  # we calculate the mean of counts
        
            err = np.sqrt(mean) / mean  # we calculate the poisson error of counts
        
            mean_level.append(mean)
        
            data_err.append(err)
        
    if mask_path != None:
        
        mask_path = Path(mask_path)
            
        mask_array = CCDData.read(mask_path, unit='adu').data
        
        for exp_time in sorted(flat_times):
        
            # we filter them by the exposure time, this will ensure that we will use only a flat of a given time
        
            flat = calibrated_flats.files_filtered(object = 'flat', exptime=exp_time,include_path=True, combined=None)
        
            if len(flat) > 1:
                raise Exception('THERE ARE AT LEAST 2 FLAT FRAMES WITH THE SAME EXPOSURE TIME, YOU NEED ONLY ONE.')
            
            flat = CCDData.read(Path(flat[0])).data
            
            flat = np.ma.masked_array(flat, mask_array)
            
            mean = np.mean(flat)  # we calculate the mean of counts
        
            err = np.sqrt(mean)/mean   # we calculate the poisson error of counts
        
            mean_level.append(mean)
        
            data_err.append(err)  
    	
    ########################################################################################################################
    
    def relative_error(y1,y2):    # for calculating the relative error of each data_point in y axis.
        y0 = abs(np.subtract(y1,y2))
        y00 = np.divide(y0,y2)
        result = y00*100          # we get the result in percertage %
        return result
        
    ########################################################################################################################
    
    x = np.array(exptimes_flats)  # converting the list into an array
    y = np.array(mean_level)
    y_err = np.array(data_err)    # we need error for our y data # poisson error
    
    ########################################################################################################################
    
    def linear(x,a,b):            # model to be fitted. Here we use a linear function
        return a*x + b
    
    ########################################################################################################################
    
    if adu_min == None:
        
        if adu_max == None:
            
            least_squares = LeastSquares(x, y, y_err,  model = linear) # we use LeastSquares function from iminuit numpy module for fitting

            m = Minuit(least_squares, a = 1, b = 0)  # starting values for a and b

            m.migrad() # finds minimum of least_squares function
            m.hesse() # accurately computes uncertainties

            correlation = np.array(m.covariance.correlation())[1,0]  # accurately computes uncertainties
            #print(v[1,0]**2)

            plt.figure(figsize = (10,8))
            plt.errorbar(x, y, y_err, fmt="o", label="data")
            plt.plot(x, linear(x, *m.values), label="fit")


            relative_err = relative_error(y, linear(x, *m.values))

            total_re = np.mean(relative_err)
            
            # display legend with some fit info
            fit_info = [f'$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(x) - m.nfit}', f'$R^2$ = {correlation**2:.3f}']
            for p, v, e in zip(m.parameters, m.values, m.errors):
                fit_info.append(f"{p} = ${v:.3f} \\pm {e:.3f}$")

            plt.title('Linearity Test using Flat Frames', fontsize=20)
            plt.ylabel('Mean Value of Counts' , fontsize = 20)    
            plt.xlabel('Exposure Time [s]', fontsize = 20)
            plt.legend(title="\n".join(fit_info), fontsize  = 15, title_fontsize=15)
            
            
            y = list(map('{:.2f}'.format,y))
            
            relative_err = list(map('{:.3f}%'.format,relative_err))
        

            print(f'YOUR CORRELATION PARAMETER IS R = {abs(correlation):.3f}, IT SHOULD BE CLOSE TO 1.')

            print('WORKING IN THE RANGE OF 0 TO 65000 ADU.')
            
            print(f'YOU DATA AFTER USING THE SELECTED RANGE OF COUNTS VALUES IS: DATA = {y} ADU.')

            print(f'THE RELATIVE ERROR OF EACH DATA POINTS ARE: RE =  {relative_err}.' )
            
            print(f'THE TOTAL RELATIVE ERROR IS: TOTAL RE = {total_re:.3f}%.')
            
            if save == True:
                
                 plt.savefig( path1 + 'Linearity_Test.png', dpi = 600,  bbox_inches='tight');
                 
    ######################################################################################################################################             

    if adu_min != None:
        
        if adu_max != None:
            
            x = np.array(exptimes_flats)  # converting the list into an array
            y = np.array(mean_level)
            y_err = np.array(data_err)    # we need error for our y data # poisson error
  
            d = np.vstack((y, x, y_err)).T
    
            d = d[d[:,0]> adu_min]
            
            d = d[d[:,0]< adu_max]
            
            y = d[:,0]
            
            x = d[:,1]
            
            y_err = d[:,2]
    
    
            least_squares = LeastSquares(x, y, y_err,  model = linear) # we use LeastSquares function from iminuit numpy module for fitting

            m = Minuit(least_squares, a = 1, b = 0)  # starting values for a and b

            m.migrad() # finds minimum of least_squares function
            m.hesse() # accurately computes uncertainties


            print(m.covariance)
            correlation = np.array(m.covariance.correlation())[1,0]  # accurately computes uncertainties
            #print(v[1,0]**2)

            plt.figure(figsize = (10,8))
            plt.errorbar(x, y, y_err, fmt="o", label="data")
            plt.plot(x, linear(x, *m.values), label="fit")


            relative_err = relative_error(y, linear(x, *m.values))
            
            total_re = np.mean(relative_err)

            # display legend with some fit info
            fit_info = [f'Total relative error: {total_re:.3f}% ', f'$R^2$ = {correlation**2:.3f}']
            #fit_info = [f'$R^2$ = {correlation**2:.3f}']
            for p, v, e in zip(m.parameters, m.values, m.errors):
                fit_info.append(f"{p} = ${v:.3f} \\pm {e:.3f}$")

            plt.title(f'Linearity Test using Flat Frames. ADU range = [{adu_min}, {adu_max}]', fontsize=20)
            plt.ylabel('Mean Value of Counts' , fontsize = 20)    
            plt.xlabel('Exposure Time [s]', fontsize = 20)
            plt.legend(title="\n".join(fit_info), fontsize  = 15, title_fontsize=15)
            plt.xlim(min(x), max(x))

            y = list(map('{:.2f}'.format,y))
            
            relative_err = list(map('{:.3f}%'.format,relative_err))
        

            print(f'YOUR CORRELATION PARAMETER IS R = {abs(correlation):.3f}, IT SHOULD BE CLOSE TO 1.')

            print(f'WORKING IN THE RANGE FROM {adu_min } TO {adu_max} ADU.')
            
            print(f'YOU DATA AFTER USING THE SELECTED RANGE OF COUNTS VALUES IS: DATA = {y} ADU.')

            print(f'THE RELATIVE ERROR OF EACH DATA POINTS ARE: RE =  {relative_err}.' )
            
            print(f'THE TOTAL RELATIVE ERROR IS: TOTAL RE = {total_re:.3f}%.')

            
            if save == True:
                
                 plt.savefig( path1 + 'Linearity_Test.png', dpi = 600,  bbox_inches='tight');

    


def readout_noise(path_bias_folder, gain = None, warning = 'ignore'):
    
    
    # we calculate the readout noise of the camera using bias frames, the goal of doing this is because thereadout noise is a statistical error source for ccd cameras.
    
    """
    path_bias_folder: path to the folder when you place all the bias frames, recommended to take 100 frames to have a better statistc.
    
    gain: float. Gain of the camera, which is typically given in electrons/ADU. By giving this value you will get the readoutnoise from the camera and the error matrix in electron units.
     
    """
    warnings.simplefilter( warning , category=AstropyWarning)  # avoiding the warnings messages from astropy
    
    ##########################################################################################################
    
    path = Path(path_bias_folder)
    
    bias = ImageFileCollection(path)  # creating a collection of images
     
    #########################################################################################################
    
    for im in bias.hdus(combined=True):  # seek if the combined bias is already there,
        if im != None:
            #raise Exception('COMBINED BIAS IS ALREADY IN THE FOLDER!')
            warnings.warn("COMBINED BIAS ALREADY IN THE FOLDER, IT WILL BE OVERWRITTEN!")
    
    #########################################################################################################
      
    biases = bias.files_filtered(object = 'bias', include_path=True, combined =None ) # filtering the image collection in case you have something else than bias frames
    
    masterbias = ccdp.combine(biases,        # getting master bias using combine functio
                             method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e6, unit='adu')

    masterbias.meta['combined'] = True  # adding the combined header
    
    ############################################################################################################
    
    subbias_list = []   # creating an empty list to place our 2D array images after subtracting the masterbias to all the bias frames
    
    for ccd, file_name in bias.ccds( object = 'bias', combined = None,     
                                      ccd_kwargs={'unit': 'adu'} ,   # CCDData requires a unit for the image if 
                                      return_fname = True ):         # Provide the file name too.
        ccd = ccdp.subtract_bias(ccd, masterbias).data   # subtracting masterbias
        
        subbias_list.append(ccd)  # filling the empty list with the result after substracting
    
    ############################################################################################################
    
    number_bias = len(subbias_list)  # printing the number of images in the list to make sure we use all the bias iamges
    
    std_list = []
    
    mean_array = np.mean(subbias_list, axis=1)

    # Calculate absolute differences
    std_list = np.abs(np.subtract(subbias_list, mean_array[:, np.newaxis, :]))

    def mean_matrix(matrices_list):
        stacked_matrices = np.stack(matrices_list, axis=0)
        mean_matrix = np.median(stacked_matrices, axis=0)
        return mean_matrix
    
    readout_noise_matrix = mean_matrix(std_list)
    
    # print(readout_noise_matrix)
    #print(f'WE HAVE USED {number_bias} BIAS FRAMES.')
    
    
    ############################################################################################################
    
    print(f'WE HAVE USED {number_bias} BIAS FRAMES.')
    
    
    if gain != None:   # in case gain value is given to the function
        if type(gain) == float:  # checking a float have been entered as the input 
            readout_noise_matrix_e = np.median(readout_noise_matrix*gain)   # calculating the readout noise of the camera which is typically given by the manufacturer and so you can compare it
            sigma_readout_noise_matrix_e = np.std(readout_noise_matrix*gain)
            print(f'THE READ OUT OF THE CAMERA IS {readout_noise_matrix_e}[e-]')   
            print(f'THE EEROR IN THE CALCULATED READ OUT OF THE CAMERA IS {sigma_readout_noise_matrix_e}[e-]')   
        
        else:
             raise ValueError("GAIN HAS TO BE A FLOAT")  # raising error in case you have not given a float to the function in the gain argument
             
    ############################################################################################################
    print(readout_noise_matrix)
    return readout_noise_matrix  # return the error matrix which after can be used which the other source of statistical error and then added to the science images
    




    
    
def dark_current_err(path_dark_frames, path_masterbias, gain = None, warning = 'ignore'):
    
    # We calculate the statistical error associated to the dark current noise of the ccd camera by using dark frames which are calibrated by substracting the masterbias
    
    """
    path_dark_frames: path to the folder where you placed the dark frames, remember they have to have the same exposure time.
    
    path_masterbias: path to the master bias. You can use the master bias you calculated for the readout noise calculation.
    
    gain: float. Gain of the camera, which is typically given in electrons/ADU. By giving this value you will get the dark current from the camera and the error matrix in electron/seconds units. 
    """
    warnings.simplefilter( warning , category=AstropyWarning)
    
    #########################################################################################
    
    masterbias = CCDData.read(Path(path_masterbias))
    
    path = Path(path_dark_frames)
    
    darks = ImageFileCollection(path)
    
    ##########################################################################################
    
    for im in darks.hdus(combined=True, unit='adu'):  # seek if the combined bias is already there,
        if im != None:
            raise Exception('COMBINED IMAGE IN THE FOLDER, CHECK IT.')
            #warnings.warn("COMBINED IMAGE IN THE FOLDER, CHECK IT.")
            
    ##########################################################################################
    
    dark_times = set(darks.summary['exptime']) # the exposure times of your dark frames
    
    dark_times_list = list(dark_times)
    
    ###########################################################################################
    
    if len(dark_times_list) > 1:  # we are checking if there are different exposure times which should not be the case
        
        print(f'EXPOSURE TIME OF YOUR DARK FRAMES ARE {dark_times_list} SECONDS.')
        raise Exception('YOU HAVE DARK FRAMES WITH DIFFERENT EXPOSURE TIMES, PLEASE CHECK IT.')
        
    if len(dark_times_list) == 1:
    
           print(f'EXPOSURE TIME OF YOUR DARK FRAMES IS {dark_times_list[0]} SECONDS.')
        
    ##########################################################################################
    
    calibrated_darks = []
    
    for ccd, file_name in darks.ccds( object = 'dark', combined = None,     
                                         # CCDData requires a unit for the image if 
                                      return_fname = True , unit = 'adu'):         # Provide the file name too.                         
                                                                                        
        # Subtract the dark current 
        
        ccd = ccdp.subtract_bias(ccd, masterbias)
        
        calibrated_darks.append(ccd.data)
        
    calibrated_darks_len = len(calibrated_darks)
    
    #############################################################################################
    
    #std_list = []
    
    
    mean_array = np.mean(calibrated_darks, axis=1)

    # Calculate absolute differences
    std_list = np.abs(np.subtract(calibrated_darks, mean_array[:, np.newaxis, :]))

    def mean_matrix(matrices_list):
        stacked_matrices = np.stack(matrices_list, axis=0)
        mean_matrix = np.median(stacked_matrices, axis=0)
        return mean_matrix
    
    dark_current_matrix = mean_matrix(std_list)
    
    
    
    #for k in calibrated_darks:
        
    #    mean = k.mean()
        
    #    std = np.array(np.matrix([abs(i-mean) for i in k]))
        
    #    std_list.append(std)
    
    #flattened_matrices = [matrix.flatten() for matrix in std_list]
    #stacked_matrices = np.vstack(flattened_matrices)

    # Compute the element-wise median of each column
    #median_array = np.median(stacked_matrices, axis=0)

    # Reshape the median array into a matrix
    
    #median_matrix = median_array.reshape(std_list[0].shape)

    ################################################################################################

    print(f'WE HAVE USED {calibrated_darks_len} DARK FRAMES')
    
    if gain != None:
    
        if type(gain) == float:
            
            time = dark_times_list[0]
            dark_current_e_s = np.median(dark_current_matrix*gain/time)
            print(f'THE DARK CURRENT OF THE CAMERA IS {dark_current_e_s}[e-/s]')
        
        else:
             raise ValueError("GAIN HAS TO BE A FLOAT")

    ###############################################################################################3

    print(dark_current_matrix)
    return dark_current_matrix
    
    
 

def att(path_image1, path_image2):

    # Final code for calculating the attenuation light coefficient. You just need to put as input two pictures you need to compare.
    
    image1 = CCDData.read(Path(path_image1), unit = 'adu')
    image2 = CCDData.read(Path(path_image2), unit = 'adu')
    
    l1 = 1.54 # [m]
    l2 = 0.54 # [m]
    
    I1 = np.sum(image1.data.flatten())
    I2 = np.sum(image2.data.flatten())
    
    err_l1 = 0.001
    err_l2 = 0.001
    
    def I_err(data):  
        x = [i**2 for i in data]
        return np.sqrt(np.sum(x))
    
    err_I1 = I_err(image1.uncertainty.array.flatten())
    err_I2 = I_err(image2.uncertainty.array.flatten())
    
    def error_propagation(func, variables, errors):
        """
        Calculates the error propagation for a function with multiple measured variables.
        :param func: The function to propagate the error for.
        :param variables: A list of measured variables.
        :param errors: A list of errors for each measured variable.
        :return: The propagated error for the function.
        """
        # Calculate the partial derivatives of the function with respect to each variable
        partials = []
        for i in range(len(variables)):
            var = variables[i]
            err = errors[i]
            x1 = np.array(variables)
            x2 = np.array(variables)
            x2[i] += err
            partial = (func(*x2) - func(*x1)) / err
            partials.append(partial)
    
        # Calculate the propagated error using the partial derivatives and errors
        delta_f = 0
        for i in range(len(partials)):
            delta_f += (partials[i] * errors[i])**2
        delta_f = np.sqrt(delta_f)
    
        return delta_f
    
    def Lambda(x1, x2, x3, x4):
        return np.abs((x1 - x2) / (np.log(x3) - np.log(x4)))
    

    # Propagate error for the function
    delta_lambda = error_propagation(Lambda, [l1, l2, I1, I2], [err_l1, err_l2, err_I1, err_I2])

    print("Result:", Lambda(l1, l2, I1, I2), "[m]^-1")
    print("Error:", delta_lambda, "[m]^-1")


























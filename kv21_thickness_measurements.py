# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:48:20 2019

@author: ebm

kv2.1 layer thickness measurements for onl and is+os
"""

import matplotlib
matplotlib.use("TkAgg")
import argparse
import pandas as pd
from tkinter.filedialog import askdirectory
import numpy as np
import skimage.morphology as mp
import skimage.measure as sm
import skimage.transform as tf
import skimage.io as io
import skimage.filters as sf
import matplotlib.pyplot as plt
import os,csv,sys,re, glob
from shutil import rmtree as rmtree
from scipy.signal import argrelextrema as locmax
from scipy.signal import savgol_filter as smooth
import scipy.ndimage as ndi
matplotlib.style.use("ggplot")

mouse_find = re.compile('(?!(20))(?<!\d)(\d{5})(?!\d)')
eye_find = re.compile('([lr]|(right|left)\s)eye', re.IGNORECASE)
#date_find = re.compile('\d{4}-\d\d-\d\d')
#for when time mark is Xweeks where X is an integer
date_find = re.compile('\d*weeks')
#geno_find = re.compile('ccr2|arr|cx3cr1|gfp/gfp|gfp/+|het|homo|gnat2|c57bl6|lox-cre|lox', re.IGNORECASE)
#Uncomment line below for Sarah's Kv 2.1 experiments. Comment line 27
geno_find = re.compile('WT|KO|het|homo', re.IGNORECASE)
#Uncomment line below for Emily's Arrestin Experiments. Comment line 25
#geno_find = re.compile('(?<!_)Arr|DRWT')
treatment_find = re.compile('saline|ccl2|clod(ronate)?', re.IGNORECASE)
ear_find = re.compile('([rlnb]|(right|left|both|neither)\s)ear', re.IGNORECASE)

path_find = re.compile('_\d\.')

def procPath(path):
    eym = eye_find.search(path)
    eam = ear_find.search(path)
    dm = date_find.search(path)
    mm = mouse_find.search(path)
    gm = geno_find.search(path)
    tm = treatment_find.search(path)

    try: date = dm.group()
    except: date = ''
    try: mouse = mm.group()
    except: mouse = ''
    try: ear = eam.group()
    except: ear = ''
    try: eye = eym.group()
    except: eye = ''
    try: treatment = tm.group()
    except: treatment = ''
    try: geno = gm.group()
    except: geno = ''
    
    date = int(date[:-5])

    return [date, mouse.lower(), ear.lower(), eye.lower(), treatment.lower(), geno.lower()]
    
def path_checker(path, i=0):
    """
    This function checks if a file exists and adds a number as an identifier at the end.
    If many functions of the same name may be made, it will increase the identifier.
    WARNING: In the case of the OCT images, a non-unique file name does NOT indicate the image is from the same animal,
    if descriptive file names were missing, the summary file names will be quite bare and could easily have
    the same names.
    Please name things well.
    """
    if os.path.isfile(path):
        if i==0:
            i+=1
            split_path = str.rsplit(path, '.')
            new_path = path_checker(split_path[0]+"_"+str(i)+"."+split_path[1], i)
        else:
            i+=1
            #use regex
            new_path = path_checker(path_find.sub('_'+str(i)+'.', path),i)
        return new_path
    else:
        return path


def findOCT(top_path):
    print(top_path)
    info_labels = ["Age (Weeks)", "mouse", "ear", "eye", "treatment", "geno","Ymin","Ymax","PR Thickness (um)","IPL Thickness (um)","IS+OS Thickness"]
    summary_path = top_path+os.sep+"SummaryInfo"
    try:
        os.mkdir(summary_path)
    except:
        rmtree(summary_path, ignore_errors=False)    
        os.mkdir(summary_path)
    with open(summary_path+os.sep+"summary_info.csv", 'w') as f:
        csvf = csv.writer(f,delimiter=",")
        csvf.writerow(info_labels)
    for path,direc,files in os.walk(top_path):
        for each in direc:
            new_path = path+os.sep+each
            if os.path.isdir(new_path) and "processed_OCT_images" in each:
                check_path = new_path+os.sep+"Log_Amp_FFT2X"
                #This is stupid and is a bandaid to process old pvOCT images.
                i_flag = False
                if os.path.isdir(check_path):
                    pass
                elif os.path.isdir(os.path.join(new_path, "I")):
                    check_path = os.path.join(new_path, "I")
                    i_flag = True
                print("Processing: ",check_path)
                if os.path.isdir(check_path):
                    if os.path.isfile(check_path+os.sep+"flat_Bscans.tif") \
                    and os.path.isfile(check_path+os.sep+"enface.tif"):
                        #get info from directory tree for labeling
                        path_info = procPath(check_path)
                        #load images
                        enface = io.imread(check_path+os.sep+"enface.tif")
                        bscans = io.imread(check_path+os.sep+"flat_Bscans.tif")
                        if i_flag is True:
                            bscans = np.flipud(bscans)
                            enface = np.flipud(enface)
                        try:
                            #find onh, save location
                            ycord = locateONH(enface)
                            #cut bscan, save
                            cut_scans = bscan_cut(bscans, ycord)
                            io.imsave(check_path+os.sep+"cut_bscan.tif", cut_scans)
                        except: pass
                        #do arrestin analysis
                        sumFig, measurements = prThickness(enface)
                        plt.tight_layout()
                        save_path = path_checker(summary_path+os.sep+str(path_info[0])+path_info[5]+path_info[2]+path_info[1]+path_info[3]+"_summary.png")
                        plt.title(str.rsplit(str.rsplit(save_path,os.sep,1)[1],'.',1)[0])                  
                        plt.savefig(save_path)
                        plt.close()

                        info = path_info+ycord+measurements

                        with open(check_path+os.sep+"oct_image_info.csv", 'w') as f:
                            cf = csv.writer(f,delimiter=",")
                            cf.writerow(info_labels)
                            cf.writerow(info)
                        with open(summary_path+os.sep+"summary_info.csv", 'a') as fi:
                            csvf = csv.writer(fi, delimiter=",")
                            csvf.writerow(info)
                        
                    else:
                        print("There was an processed_OCT_images folder with no images")
                else:
                    print("There was an processed_OCT_images folder with no Amp_FFT2X or I folder")
                    continue

def max_area(props):
    i_max=-1
    max_area = -1
    for i, prop in enumerate(props):
        bbx = np.array(prop.bbox)
        #gets rid of any region touching the sides - fufils "exlude on side" function of imagej Analyze particles
        if np.any(bbx==0) or np.any(bbx==256):
            continue
        else:
            #find max area
            if prop.area > max_area:
                max_area = prop.area
                i_max = i
    return i_max


def locateONH(enface):
    #average across x and y to get the mean a-scan profile
    profile = np.mean(enface, axis=(1,2))
    #find the max and get the index
    max_frame_ind = np.squeeze(np.array(np.where(profile==np.max(profile))))
    #get the max frame
    max_frame = enface[max_frame_ind]
    #gaussian blur
    gmax = sf.gaussian(max_frame, sigma=2)
    #automatic thresholding using the isodata algorithm
    isodata_thresh = sf.threshold_isodata(gmax)
    #invert the binary image to the onh can be segmented
    classified_img = gmax<isodata_thresh
    #find connected regions using 8-connectivity (a square)
    labels = sm.label(classified_img,connectivity=2)
    #find the properties of those regions
    props = sm.regionprops(labels)
    #detect which is the onh
    onh_ind = max_area(props)
    #save csv with centroid props[onh_ind].centroid

    #return bounding box indicies to use for b-scan cleaning function
    if onh_ind != -1:
        ymin = props[onh_ind].bbox[0]
        ymax = props[onh_ind].bbox[2]
        return [ymin, ymax]
    elif onh_ind == -1:
        #this means no onh spot was detected
        #preferable to fix, if possible.
        return [-1,-1]

def bscan_cut(bscans,onh_ybox=[-1,-1]):
    if onh_ybox!=[-1,-1]:
        #either no coordinates were entered, or they weren't detected
        cut_scans = bscans[onh_ybox[0]:onh_ybox[1],:,:] 
    else:
        cut_scans = bscans[120:136,:,:]
    avg_scans = np.mean(cut_scans, axis=0)
    #Email from Mayank to Eric on 8-15-2015 on Mayank and Robert's formula to resize bscans to 1 px**2
    #height = y_height*0.84
    #width = 1600 px ==> now about 2000 with 40um/deg
    height = np.round(avg_scans.shape[1]*0.835)
    width = 2000
    scan = tf.resize(avg_scans, (height,width), order=3, mode="reflect")
    #drop to 8 bit, otherwise the averaging takes a very long time. Which seems strange...
    scan_s = scan.astype("uint8")
    scan_m = sf.rank.mean(scan_s,selem=mp.disk(50))
    #Worried there are issues with top hat. May smear in noise that becomes part of the center of mass computation.
    scan_b = mp.white_tophat(scan_m,selem=mp.square(500))
    m = sm.moments(scan_b, order=1)
    y = int(np.round(m[0,1]/m[0,0]))
    ymin = y-300
    ymax=y+300
    cut_scan = scan[ymin:ymax,:]
    
    return cut_scan.astype("float32")

def prThickness(enface):
    #profile = smooth(np.mean(enface, axis=(1,2)),15, 3) this is worse
    profile = np.mean(enface,axis=(1,2))
    grad = ndi.gaussian_filter(np.gradient(profile), sigma=15)
    grad_min = locmax(grad, np.less, order=35)[0] #45
    # might need to blur to remove spurious extrema
    #gprofile = ndi.gaussian_filter1d(np.mean(enface, axis=(1,2)),5)
    maxLoc = locmax(profile, np.greater, order=7)[0] 
    #minLoc = locmax(profile, np.less, order=23)[0]
    divisor = 1.5
    img_mean = np.mean(profile)
    min_threshold = img_mean/divisor
    size=25
    #minLoc = []
    minLoc = locmax(profile, np.less, order=size)[0]
    y_min = profile[minLoc]
    minLoc = minLoc[np.where(y_min>min_threshold)]
    #Sarah look here:
    #CHANGE THIS TO 3 FOR DARK ADAPTED, 2 FOR LIGHT ADAPTED
    #number_of_minima is the number of minima to use in the length measurements.
    #If the incorrect optima is being detected, might need to change this.
    number_of_minima = 2
    #if number_of_minima==2:
        #print("Detecting Length for Light Adapted Retinas")
    #elif number_of_minima==3:
        #print("Detecting Length for Dark Adapted Retinas")
    #else:
        #print("Strange number of minima is being used.")
    while len(minLoc)>number_of_minima:
        minLoc = locmax(profile, np.less, order=size)[0]
        size+=2
        y_min = profile[minLoc]
        minLoc = minLoc[np.where(y_min>min_threshold)]
        
    #decreases noise threshold for dim images os minimum can be detected
    #while len(minLoc)==0:
    #    minLoc = locmax(profile, np.less, order=size)[0]
    #    divisor-=0.05
    #    min_threshold = img_mean/divisor
    #    y_min = profile[minLoc]
    #    minLoc = minLoc[np.where(y_min>min_threshold)]

    y_min = y_min[np.where(y_min>min_threshold)]
    y = profile[maxLoc]
    
    y_grad = profile[grad_min]
    grad_min=grad_min[np.where(y_grad>min_threshold)]
    y_grad = profile[grad_min]
    max_threshold = y_grad[-2]
    
    maxLoc = maxLoc[np.where(y>max_threshold)]
    #new for kv2.1, get rid of maxima in onl that are less intense than the minima 
    #in the inl
    #maxLoc = maxLoc[np.where(y>np.max(y_min))]
    y=profile[maxLoc]
    #commented out becuase threshold is below most intense minima
    #y = y[np.where(y>max_threshold)]


    #now old with grad?
    #Check for bad maxima -- will not work with bad minima... may still need to glbur
    #if y[-1]<y_min[-1]:
    #    maxLoc = maxLoc[0:-1]
    #    y = y[0:-1]

    #Below I try to select the RPE as a reference. If there is only one bright peak, 
    #it's probably not the RPE, but should be close, so I use that.
    #actually intensity, NOT location. Badly named variable.
    #RPE and bruchs? are two most intense extrema
    #rpeLoc = np.sort(y)[-2:]
    #all other extrema are not rpe
    #otherLoc = np.sort(y)[:-2]
    #they should be very close in values
    #intRatio = rpeLoc[1]/rpeLoc[0]
    #find the locations of the points in depth
    #inty = y[np.where(np.in1d(y,rpeLoc))]
    #If these two points are rpe and bruchs, the should be more intense than all other locations
    #they should be bright in value
    #and bruchs is usually brighter, if not want the brightest peak to get the RPE
    #2017-9-7: going to use bruchs because it is the most consistently present peak.
    #if np.all(rpeLoc[0]>otherLoc) and np.all(rpeLoc[1]>otherLoc) and intRatio<1.1 and inty[0]<inty[1]:
        #ref_point = maxLoc[np.where(y==rpeLoc[0])]
   # else:
        #ref_point = maxLoc[np.where(y==rpeLoc[1])]
        
        
    #OPL
    opl = maxLoc[np.where(maxLoc<grad_min[-2])][-1]
    #elm
    elm = maxLoc[np.where(maxLoc<minLoc[0])][-1]
    #rpe
    rpe = maxLoc[1]
    #ipl
    inl = minLoc[-1]
    #nfl
    nfl = grad_min[-1]
    
    pr_thickness = (opl-elm)*0.835
    is_os_thickness = (elm-rpe)*0.835
    ipl_thickness = (nfl-inl)*0.835
    
    #pr_thickness = minLoc[-1]-ref_point[0]
    #added for mid point average in PR to adjust for using bruchs instead of rpe
    #mid_pr = pr_thickness/2
    #pr_thickness = pr_thickness*0.835#convert to microns
    #ipl_thickness = (maxLoc[-1]-minLoc[-1])*0.84 #convert to microns
    #ipl_thickness = (grad_min[-1]-minLoc[-1])*0.835#convert to microns
    #pr_scatter = np.mean(enface[int(ref_point):int(minLoc[-1]),64:193,:])
    #pr_scatter_med = np.median(enface[int(ref_point):int(minLoc[-1]),64:193,:])
    #ipl_scatter = np.mean(enface[int(minLoc[-1]):int(maxLoc[-1]),64:193,:])
    #ipl_scatter_med = np.median(enface[int(minLoc[-1]):int(maxLoc[-1]),64:193,:])
    #ipl_scatter = np.mean(enface[int(minLoc[-1]):int(grad_min[-1]),64:193,:])
    #ipl_scatter_med = np.median(enface[int(minLoc[-1]):int(grad_min[-1]),64:193,:])
    #mid_scatter = np.mean(enface[int(ref_point[0]+mid_pr),64:193,:])
    #mid_scatter = np.mean(enface[int(ref_point[0]+mid_pr),64:193,:])
    #ratio = pr_scatter/ipl_scatter
    #ratio_med = pr_scatter_med/ipl_scatter_med
    #ratio_mid = mid_scatter/ipl_scatter
    
    bBad = np.average(enface, axis=1)
    scan = tf.resize(bBad, (bBad.shape[0],bBad.shape[1]*4), order=3)
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,8))
    ax[0].scatter(maxLoc, y, c="b")
    ax[0].scatter(minLoc, y_min, c='g')
    ax[0].scatter(grad_min, y_grad, c='k')
    ax[0].plot(profile)
    ax[1].imshow(scan, cmap='gray')
    ax[1].grid(b=False)
    #show horizontal lines of detected features on the b-scan
    ax[1].plot([0,scan.shape[1]-1], [rpe, rpe], lw=1,c='b', ls="--")
    ax[1].plot([0,scan.shape[1]-1], [elm, elm], lw=1,c='b', ls="--")
    ax[1].plot([0,scan.shape[1]-1], [opl, opl], lw=1,c='b', ls="--")
    ax[1].plot([0,scan.shape[1]-1], [inl, inl], lw=1,c='b', ls="--")
    ax[1].plot([0,scan.shape[1]-1], [nfl, nfl], lw=1,c='b', ls="--")
    #show thickness measurements
    ax[1].plot([500,500],[rpe, elm], ls="--")
    ax[1].plot([600,600],[elm, opl], ls="--")
    ax[1].plot([700,700],[inl, nfl], ls="--")
    h = scan.shape[0]
    
    ax[1].text(16,h-150,"PR Thickness: {:.2f}".format(pr_thickness.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-100, "IPL Thickness: {:.2f}".format(ipl_thickness.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-50, "IS+OS Thickness: {:.2f}".format(is_os_thickness.astype("float32")), color='w', family="monospace")
    #ax[1].text(16,h-50, "Ratio Midpoint: {:,.2f}".format(ratio_mid.astype("float32")),color='w', family="monospace")
    ax[1].text(500,h-50, "Number of Minima Used: {}".format(number_of_minima),color='w', family="monospace")
    ax[0].plot([0,scan.shape[0]-1], [min_threshold, min_threshold])

    return fig,[pr_thickness,ipl_thickness,is_os_thickness]

def aggImgs(path):
    #needs to be passed summary_path: ...\SummaryInfo
    img_names = glob.glob('{}{}*.png'.format(path, os.sep))
    try:
        imgs = io.imread_collection(img_names)  
    except:
        print('No images were found in {}'.format(path))
    else:
        try:
            stack = io.collection.concatenate_images(imgs)
            new_path = '{}{}img_compendium.tif'.format(path,os.sep)
            io.imsave(new_path, stack)
        except:
            sys.exit("Did you run the OCT flattening program?")
    

def aggData(path):
    summary_path = path+os.sep+"SummaryInfo"
    aggImgs(summary_path)
    try:
        df = pd.read_excel(summary_path+os.sep+"summary_info.xlsx")
    except:
        df = pd.read_csv(summary_path+os.sep+"summary_info.csv")
    #just in case false columns need to be added if a session was missed. 
    #empty columns prevent aggregation, just use the empty string for checking.
    #also need to drop empty rows that effect the agg run separately?
    df.dropna(axis=0, how='all',inplace=True)
    #alt_df = df.fillna('')
    df['treatment'] = df['treatment'].fillna('')
    #try:
        #df['date'] = pd.to_datetime(df['date'], format='%Y/%m/%d')
    #except:
        #print("Excel is stupid. You edited the csv file in excel and saved as a csv. The date format is screwed up.\nTo fix, select all the dates and right click, select 'Format Cells..'and choose the second to last option that is m/d/y with ALL 4 digits for year. Then save as xlsx. Excel is a bad program and microsoft should feel bad.")
        #sys.exit()
    #try: 
        #fails if mouse is filled with numbers, then want to use mouse
        #other wise, if any mouse fields are empty, use ear
        #assert np.any(alt_df['mouse']=='')
        #if np.all(alt_df['ear']!='')==True:
            #proc on ear
            #earliest_dates = df.groupby("ear").date.min()
            #date_diff = [(df.iloc[i]["date"]-earliest_dates[item])/pd.Timedelta(hours=1) for i,item in enumerate(df["ear"])]
            #df["Hours After Damage"] = date_diff
       # else:
            #sys.exit('\nAll mouse numbers are missing, but there is also ear information missing. You cannot tell what data will belong to which mouse.\
                    #\nProgram Terminated.')
   # except:
        #proc on mouse
        
        #any instead of all on lie 326 in the assert above should fix
        #only issue is that this will happen if any mouse numbers are present, but doesn't check that they are ALL present
        #earliest_dates = df.groupby("mouse").date.min()
        #date_diff = [(df.iloc[i]["date"]-earliest_dates[item])/pd.Timedelta(hours=1) for i,item in enumerate(df["mouse"])]
        #df["Hours After Damage"] = date_diff

    df.set_index(["Age (Weeks)","geno","treatment", "mouse","ear"], inplace=True)
    sum_df = df.groupby(level=["Age (Weeks)","geno","treatment"]).agg(["mean","sem"])
    sum_df.to_csv(summary_path+os.sep+"agg_data.csv")


def main():
    parser = argparse.ArgumentParser(description='Select Programs to Run')
    agg = parser.add_mutually_exclusive_group(required=False)
    agg.add_argument("-agg", action="store_true",help='run data aggregation only')
    agg.add_argument("-nagg", action="store_true", help='runs everything but data aggregation')
    args = parser.parse_args()
    top = askdirectory()
    if top=='':
        sys.exit("\n\nExited: No File Path Selected\n\n")
    if args.agg:
        aggData(top)
    elif args.nagg:
        findOCT(top)
    else:
        findOCT(top)
        aggData(top)
    print("Complete")

main()

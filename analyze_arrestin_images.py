import matplotlib
matplotlib.use("TkAgg")
from tkinter.filedialog import askdirectory
import numpy as np
import skimage.morphology as mp
import skimage.measure as sm
import skimage.transform as tf
import skimage.io as io
import skimage.filters as sf
import matplotlib.pyplot as plt
import os,csv,sys,re
from scipy.signal import argrelextrema as locmax
from scipy.signal import savgol_filter as smooth
import scipy.ndimage as ndi
matplotlib.style.use("ggplot")

mouse_find = re.compile('(?!(20))(?<!\d)(\d{4})(?!\d)')
eye_find = re.compile('([lr]|(right|left)\s)eye', re.IGNORECASE)
date_find = re.compile('\d{4}-\d\d-\d\d')
geno_find = re.compile('ccr2|arr|cx3cr1|gfp/gfp|gfp/+|het|homo', re.IGNORECASE)
treatment_find = re.compile('saline|ccl2|clod(ronate)?', re.IGNORECASE)
ear_find = re.compile('([rlnb]|(right|left|both|neither)\s)ear', re.IGNORECASE)

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

    return [date, mouse, ear, eye, treatment, geno]
    
def findOCT(top_path):
    print(top_path)
    info_labels = ["date", "mouse", "ear", "eye", "treatment", "geno","Ymin","Ymax","PR Thickness (µm)","IPL Thickness (µm)","Midpoint Scatter","PR Scatter","Median PR Scatter","IPL Scatter","Median IPL Scatter","Ratio PR/IPL", "Median Scatter Ratio", "Midpoint:IPL Average Scatter"]
    summary_path = top_path+os.sep+"SummaryInfo"
    try:
        os.mkdir(summary_path)
    except: pass
    with open(summary_path+os.sep+"summary_info.csv", 'w') as f:
        csvf = csv.writer(f,delimiter=",")
        csvf.writerow(info_labels)
    for path,direc,files in os.walk(top_path):
        for each in direc:
            new_path = path+os.sep+each
            if os.path.isdir(new_path) and "processed_OCT_images" in each:
                check_path = new_path+os.sep+"Amp_FFT2X"
                print("Processing: ",check_path)
                if os.path.isdir(check_path):
                    if os.path.isfile(check_path+os.sep+"flat_Bscans.tif") \
                    and os.path.isfile(check_path+os.sep+"enface.tif"):
                        #get info from directory tree for labeling
                        path_info = procPath(check_path)
                        #load images
                        enface = io.imread(check_path+os.sep+"enface.tif")
                        bscans = io.imread(check_path+os.sep+"flat_Bscans.tif")
                        #find onh, save location
                        ycord = locateONH(enface)
                        #cut bscan, save
                        cut_scans = bscan_cut(bscans, ycord)
                        io.imsave(check_path+os.sep+"cut_bscan.tif", cut_scans)
                        #do arrestin analysis
                        sumFig, measurements = prThickness(enface)
                        plt.tight_layout()
                        plt.savefig(summary_path+os.sep+path_info[0]+path_info[2]+path_info[1]+path_info[3]+"_summary.pdf")
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
                    print("There was an processed_OCT_images folder with no Amp_FFT2X folder")
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
    #width = 1600 px
    height = np.round(avg_scans.shape[1]*0.84)
    width = 1600
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
    # might need to blur to remove spurious extrema
    #gprofile = ndi.gaussian_filter1d(np.mean(enface, axis=(1,2)),5)
    maxLoc = locmax(profile, np.greater, order=7)[0] 
    minLoc = locmax(profile, np.less, order=23)[0]
    y = profile[maxLoc]
    y_min = profile[minLoc]
    threshold = np.mean(profile)/1.25
    minLoc = minLoc[np.where(y_min>threshold)]
    y_min = y_min[np.where(y_min>threshold)]
    maxLoc = maxLoc[np.where(y>threshold)]
    y = y[np.where(y>threshold)]
    
    #Below I try to select the RPE as a reference. If there is only one bright peak, 
    #it's probably not the RPE, but should be close, so I use that.
    #actually intensity, NOT location. Badly named variable.
    #RPE and bruchs? are two most intense extrema
    rpeLoc = np.sort(y)[-2:]
    #all other extrema are not rpe
    otherLoc = np.sort(y)[:-2]
    #they should be very close in values
    intRatio = rpeLoc[1]/rpeLoc[0]
    #find the locations of the points in depth
    inty = y[np.where(np.in1d(y,rpeLoc))]
    #If these two points are rpe and bruchs, the should be more intense than all other locations
    #they should be bright in value
    #and bruchs is usually brighter, if not want the brightest peak to get the RPE
    #2017-9-7: going to use bruchs because it is the most consistently present peak.
    if np.all(rpeLoc[0]>otherLoc) and np.all(rpeLoc[1]>otherLoc) and intRatio<1.1 and inty[0]<inty[1]:
        ref_point = maxLoc[np.where(y==rpeLoc[0])]
    else:
        ref_point = maxLoc[np.where(y==rpeLoc[1])]
    
    pr_thickness = minLoc[-1]-ref_point[0]
    #added for mid point average in PR to adjust for using bruchs instead of rpe
    mid_pr = pr_thickness/2
    pr_thickness = pr_thickness*0.84 #convert to microns
    ipl_thickness = (maxLoc[-1]-minLoc[-1])*0.84 #convert to microns
    pr_scatter = np.mean(enface[int(ref_point):int(minLoc[-1]),64:193,:])
    pr_scatter_med = np.median(enface[int(ref_point):int(minLoc[-1]),64:193,:])
    ipl_scatter = np.mean(enface[int(minLoc[-1]):int(maxLoc[-1]),64:193,:])
    ipl_scatter_med = np.median(enface[int(minLoc[-1]):int(maxLoc[-1]),64:193,:])
    mid_scatter = np.mean(enface[int(ref_point[0]+mid_pr),64:193,:])
    ratio = pr_scatter/ipl_scatter
    ratio_med = pr_scatter_med/ipl_scatter_med
    ratio_mid = mid_scatter/ipl_scatter
    
    bBad = np.average(enface, axis=1)
    scan = tf.resize(bBad, (bBad.shape[0],bBad.shape[1]*4), order=3)
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,8))
    ax[0].scatter(maxLoc, y, c="b")
    ax[0].scatter(minLoc, y_min, c='g')
    ax[0].plot(profile)
    ax[1].imshow(scan, cmap='gray')
    ax[1].grid(b=False)
    for item in maxLoc:
        ax[1].plot([0,scan.shape[1]-1], [item, item], lw=1,c='b', ls="--")
    for item in minLoc:
        ax[1].plot([0,scan.shape[1]-1], [item, item], lw=1,c='g', ls="--")
    ax[1].plot([512,512],[ref_point[0], minLoc[-1]], ls="--")
    ax[1].plot([600,600],[maxLoc[-1], minLoc[-1]], ls="--")
    h = scan.shape[0]
    ax[1].plot([0,scan.shape[1]-1],[int(ref_point[0]+mid_pr),int(ref_point[0]+mid_pr)])
    ax[1].text(16,h-250,"PR Thickness: {:.2f}".format(pr_thickness.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-200, "IPL Thickness: {:.2f}".format(ipl_thickness.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-150,"PR Scatter: {:,.2f}".format(pr_scatter.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-100, "IPL Scatter: {:,.2f}".format(ipl_scatter.astype("float32")), color='w', family="monospace")
    ax[1].text(16,h-50, "Scatter PR/IPL: {:,.2f}".format(ratio.astype("float32")), color='w', family="monospace")
    ax[1].text(400,h-250, "Mid Scatter: {:,.2f}".format(mid_scatter.astype("float32")),color='w', family="monospace")
    ax[1].text(400,h-200, "Median PR Scatter: {:,.2f}".format(pr_scatter_med.astype("float32")),color='w', family="monospace")
    ax[1].text(400,h-150, "Median IPL Scatter: {:,.2f}".format(ipl_scatter_med.astype("float32")),color='w', family="monospace")
    ax[1].text(400,h-100, "Ratio Median: {:,.2f}".format(ratio_med.astype("float32")),color='w', family="monospace")
    ax[1].text(400,h-50, "Ratio Midpoint: {:,.2f}".format(ratio_mid.astype("float32")),color='w', family="monospace")
    ax[0].plot([0,scan.shape[0]-1], [np.mean(profile)/1.5, np.mean(profile)/1.5])

    return fig,[pr_thickness,ipl_thickness,mid_scatter,pr_scatter,pr_scatter_med,ipl_scatter,ipl_scatter_med,ratio,ratio_med,ratio_mid]


def main():
    top = askdirectory()
    if top=='':
        sys.exit("\n\nExited: No File Path Selected\n\n")
    findOCT(top)
    print("Complete")

main()

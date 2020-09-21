"""Function for aponeuroses location"""

"*****************************************************************************"
"**********************************PACKAGES***********************************"
"*****************************************************************************"

import cv2
import numpy as np
from skimage.transform import radon, iradon

"*****************************************************************************"
"*********************************FUNCTIONS***********************************"
"*****************************************************************************"

def twoApoLocation(I, calibV, angle1, angle2,thresh = None):
    """Function that determines the radon transform of image I and segments it
    to detect the two aponeuroses as the two largest white lines.
    It returns the inverse transform of the segmented radon transform as
    well as location of two horizontal bands containing aponeuroses.
    
    Args:
        I (array): one canal image
        thresh: threshold used for segmentation. In the radon transform, all
        pixels where value > thresh are kept, so that to keep only whiter areas
        and remove gray areas. If None, the threshold is set so that it only keeps
        the 0.8% brightest pixels
        calibV: is the calibration factor along axis 0.
        
    Returns:
        (a1, b1), (a2, b2): parameters of the line equation of each
        aponeurosis: x = a * y + b where x = coordinate along axis 0 and y =
        coordinate along axis 1. Equations parameters are relative to I shape.
        loc1 (tuple): indicates two indices (distant of 50 pixels) corresponding 
        to the lines  of I between which the upper aponeurosis is.
        loc2 (tuple): indicates two indices (distant of 50 pixels) corresponding 
        to the lines  of I between which the lower aponeurosis is.
    """

    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY) 

    #working with a square because skimage radon function works on a circle only
    if I.shape[0] != I.shape[1]:
        mini = np.min(I.shape)
        I = I[0:mini,0:mini]
    
    #Calculate radon transform
    #! I_radon reverse the lines: upper aponeurosis region in I is in the bottom of
    # I_radon. Same for deep apo
    I_radon = radon(I, circle = True)
        
    #erase lines with slope > 100° and <80°
    I_radon[:, :angle1] = np.zeros((I_radon.shape[0], angle1))
    I_radon[:, angle2:] = np.zeros((I_radon.shape[0], 180-angle2))
    
    
    #cut radon transform in halves
    R_top = np.copy(I_radon[:int(I_radon.shape[0]/2),:])
    R_inf = np.copy(I_radon[int(I_radon.shape[0]/2):,:])
    
    #erase lines that are too centered in the radon transform to be aponeuroses
    fifth_x_half = int(I_radon.shape[0]/2/5)
    quarter_x_half = int(I_radon.shape[0]/2/4)
    R_top[3*fifth_x_half:,:]=np.zeros((R_top.shape[0]-3*fifth_x_half, R_top.shape[1]))
    R_inf[:quarter_x_half,:] = np.zeros((quarter_x_half,R_inf.shape[1]))
    
    #threshold to keep whitest regions
    if thresh is None:
        thresh = np.percentile(I_radon, 99.2)
    R_top = cv2.threshold(R_top, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8)
    R_inf = cv2.threshold(R_inf, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8)

    
    #enhance separation between white regions thanks to erosion
    SE = np.uint8(np.array([[0,1,0],[1,1,1],[0,1,0]]))
    R_top = cv2.erode(src = R_top, kernel = SE)
    R_inf = cv2.erode(src = R_inf, kernel = SE)
    
    #find where are white regions and contour them
    contours_top = cv2.findContours(R_top,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[0]
    contours_inf = cv2.findContours(R_inf,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[0]
    
    #find deep aponeurosis
    if len(contours_top) <1:
        return ('error', 'error'),('error', 'error'),'error','error'
    elif len(contours_top)>1:

        if len(contours_top) == 2:
            contours_top2 = contours_top
        
        #in case there are more than 2 contours,
        # if a contour is too small compared to the others, remove it.
        elif len(contours_top) > 2:
            mean_size_top = 0
            to_remove = []
            maximums = []
            
            for cont0 in range(len(contours_top)):
                mean_size_top = mean_size_top + contours_top[cont0].shape[0]
                list_radonvalues = [I_radon[contours_top[cont0][ind,0,1],\
                                            contours_top[cont0][ind,0,0]] for ind in range(contours_top[cont0].shape[0])]
                maximums.append(max(list_radonvalues))
                
            mean_size_top = mean_size_top / len(contours_top)
            marg_size = 5 #margin = 5 points in the contour 
            
            for cont0 in range(len(contours_top)):
                if contours_top[cont0].shape[0] < mean_size_top - marg_size:
                    to_remove.append(cont0)
                #if maximum value of radon transform in the contour is too small compared to the
                #others
#                elif contours_top[cont0].shape[0] >= mean_size_top - marg_size and 
#                    if maximums[cont0] < mean_of_maxRadon - 20/100*mean_of_maxRadon:
#                        to_remove.append(cont0)
#                        print('cd on maxi radon value')


            contours_top2 = [contours_top[ind] for ind in range(len(contours_top)) if ind not in to_remove]    
            
        if len(contours_top2) == 1:
            contours_top = contours_top2
        elif len(contours_top2) >= 2:
            region1 = contours_top2[0]
            region2 = contours_top2[1]
            for cont0 in range(2,len(contours_top2)):
                if np.amin(region1[:,0,1]) < np.amin(region2[:,0,1]):
                    if np.amin(contours_top2[cont0][:,0,1]) < np.amin(region2[:,0,1]):
                        region2 = contours_top2[cont0]
                if np.amin(region1[:,0,1]) >= np.amin(region2[:,0,1]):
                    if np.amin(contours_top2[cont0][:,0,1]) < np.amin(region1[:,0,1]):
                        region1 = contours_top2[cont0]                   


            # look for the closest first two contours to line 0
            # if they are less than 5 mm far from each other, then take the 2nd one
            # if they are more than 5 mm far from each other, it means they are not
            #  skin and apo, so take the first one.
            min1 = np.amin(region1[:,0,1])
            max1 = np.amax(region1[:,0,1])
            min2 = np.amin(region2[:,0,1])
            max2 = np.amax(region2[:,0,1])
            mini_dist_mm = 4
            if min1 < min2  and  abs(min2-max1)*calibV < mini_dist_mm:
                #skin and apo ; take the farthest to line 0
                #recreate contours_top with only the right region
                contours_top = [region2]
            elif min2 < min1 and abs(min1-max2)*calibV < mini_dist_mm:
                contours_top = [region1]
            elif abs(min1-max2)*calibV >= mini_dist_mm and abs(min2-max1)*calibV >= mini_dist_mm: 
                #take the closest to line 0
                if min1 < min2:
                    contours_top = [region1]
                elif min2 < min1:
                    contours_top = [region2]
                
    #find superficial aponeurosis
    if len(contours_inf) <1:
        return ('error', 'error'), ('error', 'error'), 'error', 'error'
    elif len(contours_inf)>1:

        if len(contours_inf) == 2:
            contours_inf2 = contours_inf
            
        #in case there are more than 2 contours,
        # if a contour is too small compared to the others, remove it.
        elif len(contours_inf) > 2:
            mean_size_inf = 0
            to_remove = []
            maximums = []
            
            for cont0 in range(len(contours_inf)):
                mean_size_inf = mean_size_inf + contours_inf[cont0].shape[0]
                list_radonvalues = [I_radon[contours_inf[cont0][ind,0,1], contours_inf[cont0][ind,0,0]] for ind in range(contours_inf[cont0].shape[0])]
                maximums.append(max(list_radonvalues))
                
            mean_size_inf = mean_size_inf / len(contours_inf)
            marg_size = 5 #margin = 5 points in the contour 
            
            
            for cont0 in range(len(contours_inf)):
                if contours_inf[cont0].shape[0] < mean_size_inf - marg_size:
                    to_remove.append(cont0)
                #if maximum value of radon transform in the contour is too small compared to the
                #others
#                else:
#                    if maximums[cont0] < mean_of_maxRadon - 20/100*mean_of_maxRadon:
#                        to_remove.append(cont0)
#                        print('cd on maxi radon value')
            
            
            contours_inf2 = [contours_inf[ind] for ind in range(len(contours_inf)) if ind not in to_remove]    
            
            
        #find regions closest to last line of I
        if len(contours_inf2) ==1:
            contours_inf = contours_inf2
        elif len(contours_inf2) >= 2:
            region1 = contours_inf2[0]
            region2 = contours_inf2[1]
            for cont0 in range(2,len(contours_inf2)):
                if np.amax(region1[:,0,1]) < np.amax(region2[:,0,1]):
                    if np.amax(contours_inf2[cont0][:,0,1]) > np.amax(region1[:,0,1]):
                        region1 = contours_inf2[cont0]
                if np.amax(region1[:,0,1]) >= np.amax(region2[:,0,1]):
                    if np.amax(contours_inf2[cont0][:,0,1]) > np.amax(region2[:,0,1]):
                        region2 = contours_inf2[cont0]                   

            # look for the closest first two contours to line 0
            # if they are less than 5 mm far from each other, then take the 2nd one
            # if they are more than 5 mm far from each other, it means they are not
            #  skin and apo, so take the first one.
            min1 = np.amin(region1[:,0,1])
            max1 = np.amax(region1[:,0,1])
            min2 = np.amin(region2[:,0,1])
            max2 = np.amax(region2[:,0,1])
            mini_dist_mm = 4
            if max1 < max2  and  abs(min2-max1)*calibV < mini_dist_mm:
                #skin and apo ; take the farthest to line 0
                #recreate contours_inf with only the right region
                contours_inf = [region1]
            elif max2 < max1 and abs(min1-max2)*calibV < mini_dist_mm:
                contours_inf = [region2]
            elif abs(min1-max2)*calibV >= mini_dist_mm and abs(min2-max1)*calibV >= mini_dist_mm: 
                #take the closest to line 0
                if max1 < max2:
                    contours_inf = [region2]
                elif max2 < max1:
                    contours_inf = [region1]
        
    #! offset in R_inf when search for contours
    contours_inf[0][:,0,1] = int(I_radon.shape[0]/2) + contours_inf[0][:,0,1]
    
    
    I_radonF = np.zeros(I_radon.shape)
    center_top, radius_top = cv2.minEnclosingCircle(contours_top[0][:,0,:])
    center_inf, radius_inf = cv2.minEnclosingCircle(contours_inf[0][:,0,:])
    I_radonF[int(center_top[1]), int(center_top[0])] = 255
    I_radonF[int(center_inf[1]), int(center_inf[0])] = 255
    linearApo = (iradon(I_radonF)>0)*255.
    
    

#    cv2.imshow('radon final', I_radonF)    
#    cv2.imshow('radon2', I_radon2)
#    cv2.imshow('radon3', I_radon3)
#    cv2.imshow('radon4', I_radon4)
#    cv2.imshow('inverse transform', linearApo)
#    cv2.imshow('initial image', I)
#    cv2.waitKey(0) & 0xFF
#    cv2.destroyAllWindows()

    #Determine horizontal bands containing aponeuroses
    j=0
    while j< linearApo.shape[0] and linearApo[j,int(linearApo.shape[1]/2)]==0:
        j = j+1
    upLine = max(0, j-30)
    
    j=0
    while linearApo.shape[0]-1-j>=0 and linearApo[linearApo.shape[0]-1-j,int(linearApo.shape[1]/2)]==0:
        j=j+1
    lowLine = min(linearApo.shape[0]-1-j + 30, linearApo.shape[0]-1)
    
    loc1 = (upLine,min(upLine+60, linearApo.shape[0]-1))
    loc2 = (max(0,lowLine-60),lowLine)
    
    #Determine line equation of each aponeurosis
    #ay+b=x where x is the coordinate along axis 0 and y along axis 1
    line1 = [[u,v] for u in range(linearApo[loc1[0]:loc1[1],:].shape[0])\
        for v in range(linearApo.shape[1]) if linearApo[loc1[0]:loc1[1],:][u,v]>0]
    if line1:
        line1.sort()
        if (line1[-1][1]-line1[0][1]) != 0:
            a1 = (line1[-1][0]-line1[0][0]) / (line1[-1][1]-line1[0][1])
        else:
            a1 = 0
        b1 = -a1 * line1[0][1] + line1[0][0] + loc1[0]
    else:
        a1 = 'error'
        b1 = 'error'
    
    line2 = [[i,j] for i in range(linearApo[loc2[0]:loc2[1],:].shape[0])\
         for j in range(linearApo.shape[1]) if linearApo[loc2[0]:loc2[1],:][i,j]>0]
    if line2:
        line2.sort()
        if (line2[-1][1]-line2[0][1]) != 0:
            a2 = (line2[-1][0]-line2[0][0]) / (line2[-1][1]-line2[0][1])
        else:
            a2 = 0
        b2 = -a2 * line2[0][1] + line2[0][0] + loc2[0]
    else:
        a2 = 'error'
        b2 = 'error'

    return (a1, b1), (a2, b2), loc1, loc2





def oneApoLocation(I, calibV, angle1, angle2, thresh = None):
    """
    Function that


    Args


    Returns

    Example
    
    """

    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY) 
    
    #working with a square because skimage radon function works on a circle only
    if I.shape[0] != I.shape[1]:
        mini = np.min(I.shape)
        I = I[0:mini,0:mini]
    
    #Calculate radon transform
    I_radon = radon(I, circle = True)
     
    #erase lines with slope > 100° and <80°
    I_radon[:, :angle1] = np.zeros((I_radon.shape[0], angle1))
    I_radon[:, angle2:] = np.zeros((I_radon.shape[0], 180 - angle2))

    #threshold to keep whitest regions
    if thresh is None:
        thresh = np.percentile(I_radon, 99.2)
    I_radon2 = cv2.threshold(I_radon, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8)
    
    #enhance separation between white regions thanks to erosion
    SE = np.uint8(np.array([[0,1,0],[1,1,1],[0,1,0]]))
    I_radon3 = cv2.erode(src = I_radon2, kernel = SE)

    #find where are white regions and contour them
    contours = cv2.findContours(I_radon3, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[0]

    
    #find  aponeurosis
    if len(contours) <1:
        return ('error', 'error'), 'error'
    elif len(contours)>1:
        
        if len(contours) == 2:
            contours2 = contours
        
        #in case there are more than 2 contours,
        # if a contour is too small compared to the others, remove it.
        elif len(contours) > 2:
            mean_size = 0
            to_remove = []
            
            for cont0 in range(len(contours)):
                mean_size = mean_size + contours[cont0].shape[0]
                
            mean_size = mean_size / len(contours)
            marg_size = 5 #margin = 5 points in the contour 
            
            for cont0 in range(len(contours)):
                if contours[cont0].shape[0] < mean_size - marg_size:
                    to_remove.append(cont0)


            contours2 = [contours[ind] for ind in range(len(contours)) if ind not in to_remove]    
            
        if len(contours2) == 1:
            contours = contours2
        elif len(contours2) >= 2:
            region1 = contours2[0]
            region2 = contours2[1]
            for cont0 in range(2,len(contours2)):
                if np.amin(region1[:,0,1]) < np.amin(region2[:,0,1]):
                    if np.amin(contours2[cont0][:,0,1]) < np.amin(region2[:,0,1]):
                        region2 = contours2[cont0]
                if np.amin(region1[:,0,1]) >= np.amin(region2[:,0,1]):
                    if np.amin(contours2[cont0][:,0,1]) < np.amin(region1[:,0,1]):
                        region1 = contours2[cont0]                   


            # look for the closest first two contours to line 0
            # if they are less than 5 mm far from each other, then take the 2nd one
            # if they are more than 5 mm far from each other, it means they are not
            #  skin and apo, so take the first one.
            min1 = np.amin(region1[:,0,1])
            max1 = np.amax(region1[:,0,1])
            min2 = np.amin(region2[:,0,1])
            max2 = np.amax(region2[:,0,1])
            mini_dist_mm = 4
            if min1 < min2  and  abs(min2-max1)*calibV < mini_dist_mm:
                #skin and apo ; take the farthest to line 0
                #recreate contours with only the right region
                contours = [region1]
            elif min2 < min1 and abs(min1-max2)*calibV < mini_dist_mm:
                contours = [region2]
            elif abs(min1-max2)*calibV >= mini_dist_mm and abs(min2-max1)*calibV >= mini_dist_mm: 
                #take the closest to line 0
                if min1 < min2:
                    contours = [region2]
                elif min2 < min1:
                    contours = [region1]
                    
                    
    I_radonF = np.zeros(I_radon.shape)
    center, radius = cv2.minEnclosingCircle(contours[0][:,0,:])
    I_radonF[int(center[1]), int(center[0])] = 255
    linearApo = (iradon(I_radonF)>0)*255.

    #Determine horizontal band containing aponeurosis
    j=0
    while j< linearApo.shape[0] and linearApo[j,int(linearApo.shape[1]/2)]==0:
        j = j + 1
    upLine = max(0, j-30)
        
    loc = (upLine,min(upLine+60, linearApo.shape[0]-1))
    
    #Determine line equation (in I, not in the band 'loc')
    #ay+b=x where x is the coordinate along axis 0 and y along axis 1
    line = [[u,v] for u in range(linearApo[loc[0]:loc[1],:].shape[0])\
        for v in range(linearApo.shape[1]) if linearApo[loc[0]:loc[1],:][u,v]>0]
    if line:
        line.sort(key=lambda x:x[1])
        if (line[-1][1]-line[0][1]) != 0:
            a = (line[-1][0]-line[0][0]) / (line[-1][1]-line[0][1])
        else:
            a = 0
        b = -a * line[0][1] + line[0][0] + loc[0]
    else:
        a = 'error'
        b = 'error'
    
    return (a, b), loc
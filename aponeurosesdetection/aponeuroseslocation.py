"""Function for aponeuroses location"""

"*****************************************************************************"
"**********************************PACKAGES***********************************"
"*****************************************************************************"

import math as m
import cv2
import numpy as np
from skimage.transform import radon, iradon

"*****************************************************************************"
"*********************************FUNCTIONS***********************************"
"*****************************************************************************"

def twoApoLocation(I, thresh = None):
    """Function that determines the radon transform of image I and segments it
    to detect the two aponeuroses as the two largest white lines.
    It returns the inverse transform of the segmented radon transform as
    well as location of two horizontal bands containing aponeuroses.
    
    Args:
        I (array): one canal image
        thresh: threshold used for segmentation. In the radon transform, all
        pixels where value > thresh are kept, so that to keep only whiter areas
        and remove gray areas.
    
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
    I_radon = radon(I, circle = True)
     
    #erase lines with slope > 100° and <80°
    I_radon[:, :80] = np.zeros((I_radon.shape[0], 80))
    I_radon[:, 101:] = np.zeros((I_radon.shape[0], 79))
    
    #erase lines that are in the middle third of the radon transform
    third = int(I_radon.shape[0]/3)
    I_radon[third:2*third, :] = np.zeros((third, I_radon.shape[1]))

    #threshold to keep whitest regions
    if thresh is None:
        thresh = np.percentile(I_radon, 99.2)
    I_radon2 = cv2.threshold(I_radon, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8)
    
    #enhance separation between white regions thanks to erosion
    SE = np.uint8(np.array([[0,1,0],[1,1,1],[0,1,0]]))
    I_radon3 = cv2.erode(src = I_radon2, kernel = SE)
    
    #find where are white regions and contour them
    contours = cv2.findContours(I_radon3,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[0]
    
    #Calculate bounding rectangle perimeter for each white region
    contours_tuples = []
    for i in range(len(contours)):
        p1, p2, height, width = cv2.boundingRect(contours[i]) #rectangles identified with point (p1,p2) and vectors (l2,0), (0,l1)
        #create list with tuples (i, P_rect): i is the ID of the contour, P_rect the
        #perimeter of the rectangle that bounds contour[i]
        contours_tuples.append([i, 2.*height+2.*width])
    
    #Verify that more than 2 regions have been spotted, and sort them according to size
    #-> aponeuroses = biggest white regions
    if len(contours_tuples)<2:
        raise TypeError('Less than two aponeuroses have been located. Try a lower threshold for binarization.')
    elif len(contours_tuples)>2: 
        contours_tuples.sort(key=lambda contours: contours[1])
    
    #Keep middle point of white regions to have a single line in inverse radon transform
    I_radon4 = np.zeros(I_radon3.shape)
    for x in range(len(contours_tuples)-2, len(contours_tuples)):
        id_contour = contours_tuples[x][0]
        center, radius = cv2.minEnclosingCircle(contours[id_contour][:,0,:])
        I_radon4[int(center[1]), int(center[0])] = 255
    linearApo = (iradon(I_radon4)>0)*255.
    
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
        raise ValueError('it seems that upper aponeurosis linear approximation cannot be found')
    
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
        raise ValueError('it seems that upper aponeurosis linear approximation cannot be found')
    
    return (a1, b1), (a2, b2), loc1, loc2




def oneApoLocation(I, thresh = None):
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
    I_radon[:, :80] = np.zeros((I_radon.shape[0], 80))
    I_radon[:, 101:] = np.zeros((I_radon.shape[0], 79))

    #threshold to keep whitest regions
    if thresh is None:
        thresh = np.percentile(I_radon, 99.2)
    I_radon2 = cv2.threshold(I_radon, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8)
    
    #enhance separation between white regions thanks to erosion
    SE = np.uint8(np.array([[0,1,0],[1,1,1],[0,1,0]]))
    I_radon3 = cv2.erode(src = I_radon2, kernel = SE)

    #find where are white regions and contour them
    contours = cv2.findContours(I_radon3, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[0]
    
    #Calculate bounding rectangle perimeter for each white region
    contours_tuples = []
    for i in range(len(contours)):
        p1, p2, height, width = cv2.boundingRect(contours[i]) #rectangles identified with point (p1,p2) and vectors (l2,0), (0,l1)
        #create list with tuples (i, P_rect): i is the ID of the contour, P_rect the
        #perimeter of the rectangle that bounds contour[i]
        contours_tuples.append([i, 2.*height+2.*width])
    
    #Verify that more at least 1 region has been spotted
    #sort according to rectangle size
    #-> aponeurosis = biggest white region
    if len(contours_tuples)<1:
        raise TypeError('No aponeurosis has been located. Try a lower threshold for binarization.')
    elif len(contours_tuples)>1: 
        contours_tuples.sort(key=lambda contours: contours[1])
    
    #Keep middle point of white region to have a single line in inverse radon transform
    I_radon4 = np.zeros(I_radon3.shape)
    id_contour = contours_tuples[-1][0]
    center, radius = cv2.minEnclosingCircle(contours[id_contour][:,0,:])
    I_radon4[int(center[1]), int(center[0])] = 255
    linearApo = (iradon(I_radon4)>0)*255.

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
        raise ValueError('it seems that upper aponeurosis linear approximation cannot be found')
    
    return (a, b), loc
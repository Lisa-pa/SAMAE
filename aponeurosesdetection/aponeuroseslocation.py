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

def apoLocation(I, thresh):
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
        linearApo (array): array of same size than I, where the lines detected
        equal 1, otherwise pixels equal 0: linear approximation of aponeuroses
        (a1, b1), (a2, b2): parameters of the line equation of each
        aponeurosis: x = a * y + b where x = coordinate along axis 0 and y =
        coordinate along axis 1
        loc1 (tuple): indicates two indices (distant of 50 pixels) corresponding 
        to the lines  of I between which the upper aponeurosis is.
        loc2 (tuple): indicates two indices (distant of 50 pixels) corresponding 
        to the lines  of I between which the lower aponeurosis is.
    """

    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    
    #working with a square because radon function working on a circle only
    if I.shape[0] != I.shape[1]:
        mini = np.min(I.shape)
        I = I[0:mini,0:mini]
    
    I_radon = radon(I, circle = True)
    I_radon2 = I_radon/np.max(I_radon)*255 #spreading values between 0 and 255 to enhance white points
    I_radon3 = cv2.threshold(I_radon2, thresh, 255, cv2.THRESH_BINARY)[1].astype(np.uint8) #keeping whitest regions

    contours = cv2.findContours(I_radon3,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[0] #find white objects
    contours_tuples = [(i, contours[i][:,0,:], contours[i].size) for i in range(len(contours))]

    if len(contours_tuples)<2:
        raise TypeError('Less than two aponeuroses have been located. Try a lower threshold for binarization.')
    elif len(contours_tuples)>2: #sort according to contour's size -> aponeuroses = bigger contour's size
        contours_tuples.sort(key=lambda contours: contours[2])
    
    'Keep middle point of white objects to have a single line in inverse radon transform'
    I_radon4 = np.zeros(I_radon3.shape)
    for x in range(len(contours_tuples)-2, len(contours_tuples)):
        center, radius = cv2.minEnclosingCircle(contours_tuples[x][1])
        I_radon4[int(center[1]), int(center[0])] = I_radon3[int(center[1]), int(center[0])]
    
    linearApo = (iradon(I_radon4)>0)*255.

    'Horizontal bands containing aponeuroses'
    j=0
    while linearApo[j,int(linearApo.shape[1]/2)]==0:
        j = j+1
    upLine = max(0, j-30)
    
    j=0
    while linearApo[linearApo.shape[0]-1-j,int(linearApo.shape[1]/2)]==0:
        j=j+1
    lowLine = min(linearApo.shape[0]-1-j + 30, linearApo.shape[0])
    
    loc1 = (upLine,min(upLine+60, linearApo.shape[0]))
    loc2 = (max(0,lowLine-60),lowLine)

    #equation of each line ay+b=x where x is the coordinate along axis 0 and y along axis 1
    line1 = [[u,v] for u in range(linearApo[loc1[0]:loc1[1],:].shape[0])\
        for v in range(linearApo.shape[1]) if linearApo[loc1[0]:loc1[1],:][u,v]>0]
    if line1:
        line1.sort()
        a1 = (line1[-1][0]-line1[0][0]) / (line1[-1][1]-line1[0][1])
        b1 = -a1 * line1[0][1] + line1[0][0] + loc1[0]
    else:
        raise ValueError('it seems that upper aponeurosis linear approximation cannot be found')

    line2 = [[i,j] for i in range(linearApo[loc2[0]:loc2[1],:].shape[0])\
         for j in range(linearApo.shape[1]) if linearApo[loc2[0]:loc2[1],:][i,j]>0]
    if line2:
        line2.sort()
        a2 = (line2[-1][0]-line2[0][0]) / (line2[-1][1]-line2[0][1])
        b2 = -a2 * line2[0][1] + line2[0][0] + loc2[0]
    else:
        raise ValueError('it seems that upper aponeurosis linear approximation cannot be found')


    return linearApo, (a1, b1), (a2, b2), loc1, loc2
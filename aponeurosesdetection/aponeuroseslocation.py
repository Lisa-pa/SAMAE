"""Main module."""

"*****************************************************************************"
"**********************************PACKAGES***********************************"
"*****************************************************************************"

import math as m
import cv2
import numpy as np
from skimage.transform import radon, iradon
import tkinter as tk
from tkinter.filedialog import askopenfilename

"*****************************************************************************"
"*********************************FUNCTIONS***********************************"
"*****************************************************************************"

def d2_gaussianMasks(s):
    """Implements Gaussian's second derivatives masks for a given standard
    deviation s.
    The masks' size is set as 6*s rounded to the superior integer (+1 if even).

    Args:
        s (double) : standard deviation

    Returns:
        mask_xx, mask_xy, mask_yy (array_like): masks obtained from derivation of a gaussian,
        according to x and x, to x and y, to y and y respectively.

    """
    
    'Calculation of the size of the masks'
    size = m.ceil(s*6)
    if size%2==0: #if size if even, make it odd
        size+=1
        
    half =  m.floor(size/2)
    mask_xx = np.zeros((size,size))
    mask_xy = np.zeros((size,size))
    mask_yy = np.zeros((size,size))

    'Second derivatives expressions'
    for x in range (-half, half +1):
        for y in range (-half, half +1):

            mask_xx[x + half][y + half] = (x**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
                                                
            mask_xy[x + half][y + half] = x*y/(2.*m.pi*s**6)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
                                                
            mask_yy[x + half][y + half] = (y**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
    return mask_xx, mask_xy, mask_yy

'-----------------------------------------------------------------------------'

def MVEF_2D(I, scales, thresholds):
    """Multiscale Vessel Enhancement Method for 2D images - based on Frangi's,
    Rana's, and Jalborg's works. This method searches for geometrical 
    structures which can be regarded as tubular.
    
    Args:
        I (2D array):       I is a grayscale image (otherwise it is 
                            converted to grayscale)
        thresholds (list):  thresholds is a list of two thresholds that 
                            control the sensitivity of the line filter 
                            to the measures Fr and R.  
        scales (list):      scales is a list of lengths that correspond to 
                            the diameter of the tube-like to find
    
    Outputs:
        I2 (2D array):  image I filtered by the multiscale vessel enhancement
                        filter
    References:
        Based on [ Automatic detection of skeletal muscle architecture 
        features, Frida Elen Jalborg, Master’s Thesis Spring 2016 ]
                [Frangi, 1998, Multiscale Vessel Enhancement filtering]
                [ Rana et al., 2009, Automated tracking of muscle fascicle
                orientation in B-mode ultrasound images]
                [ Rana et al., 2011, In-vivo determination of 3D muscle
                architecture of human muscle using free hand ultrasound]

    """
    if len(I.shape)>2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    vesselness = np.zeros((I.shape[0], I.shape[1], len(scales)))
    I2 = np.zeros((I.shape[0], I.shape[1]))
    
    for sc in range(len(scales)):    

        'Calculation of the Hessian matrix coefficients'
        mask_xx = cv2.flip(d2_gaussianMasks(scales[sc])[0],flipCode = -1)
        mask_xy = cv2.flip(d2_gaussianMasks(scales[sc])[1],flipCode = -1)
        mask_yy = cv2.flip(d2_gaussianMasks(scales[sc])[2],flipCode = -1)
        
        hessian_xx = cv2.filter2D(I, -1, mask_xx, anchor = (-1,-1))  
        hessian_xy = cv2.filter2D(I, -1, mask_xy, anchor = (-1,-1))
        hessian_yy = cv2.filter2D(I, -1, mask_yy, anchor = (-1,-1))

        frobeniusnorm =  np.zeros(I.shape)
        for u in range(frobeniusnorm.shape[0]):
            for v in range(frobeniusnorm.shape[1]):
                frobeniusnorm[u,v] = np.sqrt(hessian_xx[u,v]**2+\
                    hessian_yy[u,v]**2+2*hessian_xy[u,v]**2)
        b = thresholds[0]
        c = thresholds[1]
        if b == 0:
            raise ValueError('first element of thresholds cannot be null')
        if c == 0:
            c = 1 / 2 * np.amax(frobeniusnorm)
        print(b,c)
        for i in range(I.shape[0]):
            for j in range(I.shape[1]):
           
                'For each pixel, find and order the eigenvalues of the '
                'Hessian matrix'
                H = np.array([[hessian_xx[i,j], hessian_xy[i,j]],\
                              [hessian_xy[i,j], hessian_yy[i,j]]])
                eigvals, eigvects = np.linalg.eig(H)

                #reordering eigenvalues in increasing order if needed:
                if abs(eigvals[0])>abs(eigvals[1]):
                    eigvals[0], eigvals[1] = eigvals[1], eigvals[0]
                    eigvects[:,0],eigvects[:,1] = eigvects[:,1], eigvects[:,0]
                
                'Calculation of vesselness: looking for a POSITIVE highest'
                'eigenvalue <=> looking for dark tube-like structures; looking'
                'for a NEGATIVE highest eigen value <=> looking for bright'
                'tube-like structures'
                #bright tube-like structures search did not work so the 
                #temporary solution is to inverse the ultrasound image 
                #and search for dark tube-like structures
                if eigvals[1]<=0:
                    vesselness[i,j,sc] = 0

                if eigvals[1]>0:
                    
                    'ratio - for second order ellipsoid'
                    R = eigvals[0]/eigvals[1]

                    'Frobenius norm - for second order structureness'
                    Fr = m.sqrt(eigvals[0]*eigvals[0] + eigvals[1]*eigvals[1])

                    vesselness[i,j,sc]=scales[sc]*m.exp(-R*R/(2.*b*b))*(1.-m.exp(-Fr*Fr/(2*c*c)))

    'Keep the highest value of vesselness across all scales'
    for ind1 in range(I2.shape[0]):
        for ind2 in range(I2.shape[1]):
            I2[ind1,ind2] = np.max(vesselness[ind1,ind2,:])

    return I2
'-----------------------------------------------------------------------------'
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

# "*****************************************************************************"
# "*********************************TEST****************************************"
# "*****************************************************************************"
# image = cv2.imread('C:/Users/Lisa Paillard/Desktop/AponeurosesDetection/aponeurosesdetection/data/skmuscle.jpg', -1)
# imageG = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

# aponeuroses_Linear, equ1, equ2, apon1, apon2 = apoLocation(image, 220.)
# ypoints = np.arange(0,image.shape[1]-1,1)
# xpoints = np.int32((ypoints-equ1[1])/equ1[0])
# for ind in range(ypoints.shape[0]):
#     if xpoints[ind]<image.shape[0] and xpoints[ind]>=0:
#         image[xpoints[ind],ypoints[ind],:] = [0,127,255]
# #color in orange to see what has been spotted
# for x in range(aponeuroses_Linear.shape[0]):
#     for y in range(aponeuroses_Linear.shape[1]):
#         if aponeuroses_Linear[x,y] >0:
#             image[x,y,:] = [0,127,255]

# #cv2.imwrite('Radon_cropped_Kevin_jamon_20181002_153734_image.jpg', image);
# cv2.imshow('image ini',image)
# cv2.imshow('Trasnformed I',aponeuroses_Linear)
# cv2.imshow('apo1', image[apon1[0]:apon1[1],:])
# cv2.imshow('apo2', image[apon2[0]:apon2[1],:])
# cv2.waitKey(0) & 0xFF
# cv2.destroyAllWindows()


# "*****************************************************************************"
# "*****************************************************************************"
# "*****************************************************************************"
"""Functions for fascicle detection - FaDe"""

import math as m
import cv2
import numpy as np
#------------------------------------------------------------------------------#
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
    from skimage.feature import hessian_matrix, hessian_matrix_eigvals
    
    if len(I.shape)>2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)


    vesselness = np.zeros((I.shape[0], I.shape[1], len(scales)))
    I2 = np.zeros((I.shape[0], I.shape[1]))
    b = thresholds[0]
    c = thresholds[1]
    if b == 0:
        raise ValueError('first element of thresholds cannot be null')
        
    for sc in range(len(scales)): 
        H = hessian_matrix(image = I, sigma = scales[sc], order = 'rc')
        eigvals = hessian_matrix_eigvals(H) #2 eigenvalues in decreasing order; shape of ei = (2,I.shape[0], I.shape[1])
        frobeniusnorm =  np.sqrt(H[0]**2+ H[1]**2 + 2 * H[2]**2)
        if c == 0:
            c = 1 / 2 * np.amax(frobeniusnorm)

        for i in range(I.shape[0]):
            for j in range(I.shape[1]):

                'Calculation of vesselness: looking for a POSITIVE highest'
                'eigenvalue <=> looking for dark tube-like structures; looking'
                'for a NEGATIVE highest eigen value <=> looking for bright'
                'tube-like structures'
                #bright tube-like structures search did not work so the 
                #temporary solution is to inverse the ultrasound image 
                #and search for dark tube-like structures
                
                #if eigvals[0,i,j]<=0:
                    #vesselness[i,j,sc] = 0
                    #not necessary since vesselness initialized with zeros

                if eigvals[0,i,j]>0:
                    
                    'ratio - for second order ellipsoid'
                    R = eigvals[1,i,j]/eigvals[0,i,j]

                    'Frobenius norm - for second order structureness'
                    Fr = m.sqrt(eigvals[0,i,j]*eigvals[0,i,j] + eigvals[1,i,j]*eigvals[1,i,j])

                    vesselness[i,j,sc]=scales[sc]*m.exp(-R*R/(2.*b*b))*(1.-m.exp(-Fr*Fr/(2*c*c)))

    'Keep the highest value of vesselness across all scales'
    for ind1 in range(I2.shape[0]):
        for ind2 in range(I2.shape[1]):
            I2[ind1,ind2] = np.max(vesselness[ind1,ind2,:])
    
    return I2
'-----------------------------------------------------------------------------'
"""
from skimage.transform import radon, iradon
def normalizedRadon(I):


    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    
    #working with a square because radon function working on a circle only
    if I.shape[0] > I.shape[1]:
        I = I[int(I.shape[0]/2 - I.shape[1]/2):int(I.shape[0]/2 + I.shape[1]/2), :]
    else:
        I = I[:, int(I.shape[1]/2 - I.shape[0]/2):int(I.shape[1]/2 + I.shape[0]/2)]

    cv2.imshow('Part of the image analyzed', I)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()

    R = radon(I, circle = True)
    cv2.imshow('Radon Transform', R/np.max(R))
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    X = np.arange(0,I.shape[0],1)
    Y = np.arange(0,I.shape[1],1)

    for rho in range(R.shape[0]):
        for theta in range(R.shape[1]):
            line = rho - m.cos(theta)*X + m.sin(theta)*Y
            #line = (np.abs(rho - m.cos(theta)*X + m.sin(theta)*Y) <0.6)*1. + (np.abs(rho - m.cos(theta)*X + m.sin(theta)*Y) >= 0.6)*0.
            dirac = 1./(m.pi*(1.*1.+line*line))        

            L = np.sum(dirac)
            if L != 0:
                R[rho, theta] = R[rho, theta]/L

    cv2.imshow('Radon Transform', R/np.max(R))
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    R = cv2.threshold(R, 254., 255, cv2.THRESH_BINARY)[1]

    return R,I
"""
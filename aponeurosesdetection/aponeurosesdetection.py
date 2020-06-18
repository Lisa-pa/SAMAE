"""Main module."""

"*****************************************************************************"
"**********************************PACKAGES***********************************"
"*****************************************************************************"

import math as m;
import cv2;
import numpy as np;

"*****************************************************************************"
"*********************************FUNCTIONS***********************************"
"*****************************************************************************"

def d2_gaussianMasks(s):
    """Implements Gaussian's second derivatives masks for a given scale s, 
    which corresponds to the standard deviation"""
    
    'Calculation of the size of the masks'
    size = m.ceil(s*6);
    if size%2==0: #if size if peer, make it odd
        size+=1;
        
    half =  m.floor(size/2);
    mask_xx = np.zeros((size,size));
    mask_xy = np.zeros((size,size));
    mask_yy = np.zeros((size,size));

    'Second derivatives expressions'
    for x in range (-half, half +1):
        for y in range (-half, half +1):

            mask_xx[x + half][y + half] = (x**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2));
                                                
            mask_xy[x + half][y + half] = x*y/(2.*m.pi*s**6)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2));
                                                
            mask_yy[x + half][y + half] = (y**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2));
    return mask_xx, mask_xy, mask_yy

'-----------------------------------------------------------------------------'

def MVEF_2D(I, scales, thresholds):
    """Multiscale Vessel Enhancement Method for 2D images - based on Frangi's,
    Rana's, and Jalborg's works.
    I is a grayscale image. thresholds is a list of two thresholds which 
    control the sensitivity of the line filter to the measures Fr and R.  
    scales is a list of lengths that correspond to the diameter of the 
    tube-like structure the operator wants to find"""
    
    vesselness = np.zeros((I.shape[0], I.shape[1], len(scales)));
    I2 = np.zeros((I.shape[0], I.shape[1]));  
    
    for sc in range(len(scales)):    

        'Calculation of the Hessian matrix coefficients'
        mask_xx = cv2.flip(d2_gaussianMasks(scales[sc])[0],flipCode = -1);
        mask_xy = cv2.flip(d2_gaussianMasks(scales[sc])[1],flipCode = -1);
        mask_yy = cv2.flip(d2_gaussianMasks(scales[sc])[2],flipCode = -1);
        
        hessian_xx = cv2.filter2D(I, -1, mask_xx, anchor = (-1,-1));    
        hessian_xy = cv2.filter2D(I, -1, mask_xy, anchor = (-1,-1));
        hessian_yy = cv2.filter2D(I, -1, mask_yy, anchor = (-1,-1));
        
#        cv2.imshow('maskx',hessian_xx);
#        cv2.imshow('maskxy',hessian_xy);
#        cv2.imshow('masky',hessian_yy);
#        cv2.waitKey(0) & 0xFF;
#        cv2.destroyAllWindows();
        
        for i in range(I.shape[0]):
            for j in range(I.shape[1]):
           
                'For each pixel, find and order the eigenvalues of the '
                'Hessian matrix'
                H = np.array([[hessian_xx[i,j], hessian_xy[i,j]],\
                              [hessian_xy[i,j], hessian_yy[i,j]]]);
                eigvals, eigvects = np.linalg.eig(H);
                #reordering eigenvalues in increasing order if needed:
                if abs(eigvals[0])>abs(eigvals[1]):
                    eigvals[0], eigvals[1] = eigvals[1], eigvals[0];
                    eigvects[:,0],eigvects[:,1] = eigvects[:,1], eigvects[:,0];
                
                'Calculation of vesselness: looking for a POSITIVE highest'
                'eigenvalue <=> looking for dark tube-like structures; looking'
                'for a NEGATIVE highest eigen value <=> looking for bright'
                'tube-like structures'
                #bright tube-like structures search did not work so the 
                #temporary solution is to inverse the ultrasound image 
                #and search for dark tube-like structures
                if eigvals[1]<=0:
                    vesselness[i,j,sc] = 0;

                if eigvals[1]>0:
                    
                    'ratio - for second order ellipsoid'
                    R = eigvals[0]/eigvals[1]; 

                    'Frobenius norm - for second order structureness'
                    Fr = m.sqrt(eigvals[0]*eigvals[0] + eigvals[1]*eigvals[1]);

                    b = thresholds[0];
                    c = thresholds[1];
                    vesselness[i,j,sc]=scales[sc]*m.exp(-R*R/(2.*b*b))*(1.-m.exp(-Fr*Fr/(2*c*c)));

    'Keep the highest value of vesselness across all scales'
    for ind1 in range(I2.shape[0]):
        for ind2 in range(I2.shape[1]):
            I2[ind1,ind2] = np.max(vesselness[ind1,ind2,:]);

    return I2, vesselness, hessian_xy, hessian_xx, hessian_yy
'-----------------------------------------------------------------------------'

"*****************************************************************************"
"*********************************PROGRAM*************************************"
"*****************************************************************************"

'Opening of the image and color transform.'
RGB_image = cv2.imread('cropped0.jpg', -1);
gray_I = cv2.cvtColor(RGB_image, cv2.COLOR_RGB2GRAY);
gray_I2 = 255 - gray_I; #inverse black and white
#gray_I2 = cv2.equalizeHist(gray_I);

'Multiscale Vessel Enhancement Method'
seg, vesselness,hessXY, hessXX, hessYY = MVEF_2D(gray_I2, [4.], [0.5, 0.5]);

'Rescale the pixels values if max. value is too low'
maxi = np.max(seg);
if maxi <=0.5 and maxi>0:
    seg = (seg/maxi)*255;

'Visualization'
cv2.imshow('MVEF',seg);
cv2.imshow('Initial image', gray_I);
cv2.imshow('Inverse image', gray_I2);
cv2.waitKey(0) & 0xFF;
cv2.destroyAllWindows();

"*****************************************************************************"
"*****************************************************************************"
"*****************************************************************************"
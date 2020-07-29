import cv2
import numpy as np

def autocropping(I, threshCmin, threshCmax, threshRmin, threshRmax, calibV = 0):
    """ Cropping of raw ultrasound image I to get region of interest: Removal 
    of lateral, top and bottom strips by thresholding the mean pixel value 
    of raws and columns.
    If calibV is not zero, then 2 mm are removed from the top of I (this 
    program assumes that it corresponds to the skin location).
    
    Args:
        I (array): grayscale or RGB image (the function will convert it to 
        grayscale image if it has not yet been done)
        threshCmin (double): between 0. and 255., minimum threshold value 
        for the lateral strips
        threshCmax (double): between 0. and 255., maximum threshold value 
        for the lateral strips
        threshRmin (double): between 0. and 255., minimum threshold value 
        for the bottom and top strips
        threshRmax (double): between 0. and 255., maximum threshold value 
        for the bottom and top strips
        calibV (double): vertical calibration factor. Default value is 0.
        If calibV is different from 0, then the equivalent of 2mm is
        cropped at the top of I to remove skin from the image.

    Returns:
        I2 (array): same dimension as I, containing only the region of interest

    References:
        Based on [ Automatic detection of skeletal muscle architecture 
        features, Frida Elen Jalborg, Master’s Thesis Spring 2016 ]

    Example:
        > RGB_image = cv2.imread('AponeurosesDetection/aponeurosesdetection/data/simple_echo.jpg', -1)
          #change this example: should use the function simpleimg() to load the image simple_echo, 
          #instead of using cv2.imread with a path
          cropped_img = autocropping(RGB_image, 10., 15., 12., 25.)
          cv2.imshow('image',cropped_img)
          cv2.waitKey(0) & 0xFF
          cv2.destroyAllWindows()
"""
    #Check if I is a 1-canal image
    if len(I.shape) == 3:
        grayI = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    elif len(I.shape) == 2:
        grayI = I

    RightColumnsBelow = [0]
    LeftColumnsBelow = [0]
    UpRawsBelow = [0]
    BottomRawsBelow = [0]
    
    #vertical cropping
    while threshCmin < threshCmax and \
    (len(LeftColumnsBelow)==1 or len(RightColumnsBelow)==1):
        
        for col in range(int(grayI.shape[1]/2)):
            if np.mean(grayI[:,col])<=threshCmin:
                RightColumnsBelow.append(col)
            if np.mean(grayI[:,grayI.shape[1]-1-col])<=threshCmin:
                LeftColumnsBelow.append(grayI.shape[1]-1-col)
        threshCmin+=1.
        
    #horizontal cropping
    while threshRmin < threshRmax and \
    (len(UpRawsBelow)==1 or len(BottomRawsBelow)==1):
        
        for lig in range(int(grayI.shape[0]/2)):
            if np.mean(grayI[lig,:])<=threshRmin:
                UpRawsBelow.append(lig)
            if np.mean(grayI[grayI.shape[0]-1-lig,:])<=threshRmin:
                BottomRawsBelow.append(grayI.shape[0]-1-lig)
        threshRmin+=2.  
    
    if len(I.shape) == 3:
        I2 = I[UpRawsBelow[-1]:BottomRawsBelow[-1],RightColumnsBelow[-1]:LeftColumnsBelow[-1],:]
    elif len(I.shape) == 2:
        I2 = I[UpRawsBelow[-1]:BottomRawsBelow[-1],RightColumnsBelow[-1]:LeftColumnsBelow[-1]]
    
    if calibV != 0:
        twomm = int(2 / calibV)
        I2 = I2[twomm:,:]
        
    return I2
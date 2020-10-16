

def preprocessingApo(I, typeI, mode, margin, sizeContrast):
    """This function aims at enhancing the image I -and particularly aponeuroses-
    before it is given to the active contour function.
    It applies a median filter to I, then enhances the contrast thanks to
    a chosen mode (global contrast enhancement or local contrast enhancement -
    mean, median or midgrey local contrast enhancement). Finally the image
    undergoes a morphological opening to remove small branches and irregularities
    from the aponeurosis. The structuring element SE depends on the type of
    image (simple or panoramic)
            
        Args:
                I (np-array) : 1 canal or three-canal image. If the image has
                three canals, it is converted to grayscale image.
                mode (string): chosen contrast method
                        'global': pixels whose value is less than 120 are set
                        to zero, while pixels whose value is above 120 are kept.
                        These latter are redistributed between 0 and 255.
                        'localmedian': median is calculated for each pixel
                        in a neighborhood of sizeContrast x sizeContrast. If 
                        pixel > median - margin, then the pixel's value is 
                        unchanged. Otherwise it is set to zero.
                        'localmean': mean is calculated for each pixel
                        in a neighborhood of sizeContrast x sizeContrast. If 
                        pixel > mean - margin, then the pixel's value is 
                        unchanged. Otherwise it is set to zero.
                        'localmidgrey': maximum and minimum of a 
                        sizeContrast x sizeContrast neighborhood of each pixel
                        are found. If pixel > (max+min)/2 - margin, then the 
                        pixel's value is unchanged. Otherwise it is set to zero.
                margin (int) : constant used in contrast enhancement (see section above)
                sizeContrast (int): odd integer that caracterizes neighborhood
                size in local contrast enhancement. It should be big
                enough to have a representative neighborhood of the pixel. We
                recommend sizeContrast = 41.

        Output:
                np-array of same size than I
    """
    import numpy as np
    import cv2
    import scipy.ndimage as scindi
    
    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
        
    #-----denoising-----#
    #Non local means : denoising - cv2 function is not working#
    #I2 = cv2.fastNlMeansDenoising(I, h = 3, templateWindowSize = 7, searchWindowSize = 21, normType = cv2.NORM_L2)
    
    #-----median filter-----#
    I2 = cv2.medianBlur(I, 5)

    #-----contrast-----#
    if mode == 'global':
        I2 = (I2<120)*0. + ((I2>=120)*I2 - np.min(I2)) / (np.max(I2)-np.min(I2))
        I2 = (I2>255)*255+(I2<0)*0 + ((I2>0) * (I2<255))*I2        
    
    elif mode == 'localmedian':
        Im = cv2.medianBlur(I2, sizeContrast)
        for x in range (I2.shape[0]):
            for y in range (I2.shape[1]):
                if (int(I2[x,y])-int(Im[x,y])+int(margin))<0:
                    I2[x,y] = 0
                else:
                    I2[x,y] = int(I2[x,y])
                    
    elif mode == 'localmidgrey':
        Imax = scindi.maximum_filter(I2, size = (sizeContrast, sizeContrast), mode = 'reflect', origin = 0)
        Imin = scindi.minimum_filter(I2, size = (sizeContrast, sizeContrast), mode = 'reflect', origin = 0)
        Im = (Imax + Imin)/2
        for x in range (I2.shape[0]):
            for y in range (I2.shape[1]):
                if (int(I2[x,y])-int(Im[x,y])+int(margin))<0:
                    I2[x,y] = 0
                else:
                    I2[x,y] = int(I2[x,y])
                    
    elif mode == 'localmean':
        filt = np.ones((sizeContrast, sizeContrast))*1/(sizeContrast**2)
        Im = cv2.filter2D(I2, -1, filt, (-1,-1), borderType = cv2.BORDER_REFLECT)
        
        for x in range (I2.shape[0]):
            for y in range (I2.shape[1]):
                if (int(I2[x,y])-int(Im[x,y])+int(margin))<0:
                    I2[x,y] = 0
                else:
                    I2[x,y] = int(I2[x,y])
                    
    #-----morphological operations: opening-----#
    if typeI == 'simple':
        SE = np.uint8(np.array([[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]]))
        #1) Erosion#
        I2 = cv2.erode(src=I2, kernel = SE,anchor=(-1,-1), iterations= 3, borderType = cv2.BORDER_REPLICATE)
        #2) dilatation#
        I2 = cv2.dilate(src = I2, kernel = SE,anchor = (-1,-1), iterations=3, borderType= cv2.BORDER_REPLICATE)
       
    elif typeI == 'panoramic':
        SE = np.uint8(np.array([[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]))
        #SE = np.uint8(np.array([[0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.],[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],[0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.]]))

        #1) Erosion#
        I2 = cv2.erode(src=I2, kernel = SE,anchor=(-1,-1), iterations=3, borderType = cv2.BORDER_REPLICATE)
        #2) dilatation#
        I2 = cv2.dilate(src = I2, kernel = SE,anchor = (-1,-1), iterations=3, borderType= cv2.BORDER_REPLICATE)

    return I2
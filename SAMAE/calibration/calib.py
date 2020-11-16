import cv2
import numpy as np

def autoCalibration(I):
    """Returns horizontal and vertical factors by which every distance in
    pixels should be multiplied in order to obtain the equivalent distance in
    millimeters. This program assumes that the scale presents clear axis ticks and
    that the distance between two biggest ticks is equal to 10 mm.
    It also assumes that both horizontal and vertical scales are present in the
    up right quarter of image I.
    
    Args:
        I (array): one canal image. If I is a RGB image, it is transformed
        to a grayscale image.

    Returns:
        calibFactorX (double) and calibFactorY (double) are respectively the
        horizontal and vertical calibration factors
    """
    
    #Check if I is a 1-canal image
    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    length, width = I.shape[0], I.shape[1]
    
    #Cropping with empirical percentages and binarization of the selection
    # !!! EMPIRICAL
    TCP = 0.1 #Top cropping percentage - #empirical percentage
    LCP = 0.5 #Left cropping percentage
    BCP = 0.65 #Bottom cropping percentage
    RCP = 0.1 #Right cropping percentage
    Scale_image = I[int(TCP * length):length-int(BCP * length),\
                    int(LCP * width):width-int(RCP * width)]  

    Binar_I = cv2.threshold(Scale_image, 220., 255, cv2.THRESH_BINARY)[1]                


        
    #Selection of the biggest axis ticks: contours of white objects are found as
    #well as minimal rectangles encapsulating each object. Conditions on the
    #size of these contours/bounding rectangles enable the removal of objects
    #that are not the biggest ticks
    contours = cv2.findContours(Binar_I, cv2.RETR_EXTERNAL, \
                                           cv2.CHAIN_APPROX_NONE)[0]
    
    contours_size = [contours[i].size for i in range (len(contours))]
    BoundingRectangles = []
    
    for i in range(len(contours)):
        if contours_size[i]<=1.7*np.mean(contours_size): #condition to stop considering the objects corresponding to figures
            p1, p2, l1, l2 = cv2.boundingRect(contours[i]) #rectangles identified with point (p1,p2) and vectors (l2,0), (0,l1)
            BoundingRectangles.append([i, (p1,p2,l1,l2), 2.*l1+2.*l2])
    MeanPerim = np.mean([BoundingRectangles[i][2] for i in range(len(BoundingRectangles))])
    Dashes = [BoundingRectangles[i] for i in range(len(BoundingRectangles)) if BoundingRectangles[i][2]>MeanPerim] #removal of points and small dashes

    
    #Calculation of the minimal distances between two horizontal ticks and
    #two vertical ticks
    #browse all detected axis ticks
    horiz = 10000000.
    vertic = 10000000.
    for i in range (0, len(Dashes)-1):
        ref_Dash = Dashes[i][1]
        for j in range(i+1,len(Dashes)):
            
            if len(set(list(range(Dashes[j][1][0],Dashes[j][1][0]+Dashes[j][1][2])))\
                   .intersection(list(range(ref_Dash[0],ref_Dash[0]+ref_Dash[2]))))>2:
                
                h = abs(ref_Dash[1]+ref_Dash[3]-Dashes[j][1][1]-Dashes[j][1][3])
                if h<vertic:
                    vertic = h
                                
            if len(set(list(range(Dashes[j][1][1],Dashes[j][1][1]+Dashes[j][1][3])))\
                   .intersection(list(range(ref_Dash[1],ref_Dash[1]+ref_Dash[3]))))>2:
                
                h = abs(ref_Dash[0]-Dashes[j][1][0])
                                            
                if h<horiz:
                    horiz = h             

    #Factors to convert distance in pixels into distance in millimeters
    if horiz == 10000000. or horiz == 0:
        calibFactorX = None
    else:
        calibFactorX = 10./horiz
        
    if vertic == 10000000. or vertic == 0:
        calibFactorY = None
    else:
        calibFactorY = 10./vertic
        
    ''' visual check  
    for d in range(len(Dashes)):
        p1 = Dashes[d][1][0]
        p2 = Dashes[d][1][1]
        l1 = Dashes[d][1][2]
        l2 = Dashes[d][1][3]
        for l in range(p1,p1+l1+1):
            Binar_I[p2,l] = 150
            Binar_I[p2+l2,l] = 150
        for c in range(p2,p2+l2+1):
            Binar_I[c,p1] = 150
            Binar_I[c,p1+l1] = 150            
    cv2.imshow('Binary image', Binar_I)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()   
    '''
    return calibFactorX, calibFactorY
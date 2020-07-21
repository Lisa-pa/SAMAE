import cv2
import numpy as np

def autoCalibration(I):
    """Returns horizontal and vertical factors by which every distance in
    pixels should be multiplied in order to obtain the equivalent distance in
    millimeters. This program assumes that the scale presents clear dashes and
    that the distance between two biggest dashes is equal to 10 mm.
    
    Args:
        I (array): one canal image. If I is a RGB image, it is transformed
        to a grayscale image.

    Returns:
        calibFactorX (double) and calibFactorY (double) are respectively the
        horizontal and vertical calibration factors

    Example:
        > import tkinter as tk
          from tkinter.filedialog import askopenfilename
          root = tk.Tk()
          root.withdraw()
          filename = askopenfilename(initialdir = "/",title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
          RGB_image = cv2.imread('filename', -1)
          calibX, calibY = autoCalibration(RGB_image)
          print(calibX, calibY)
          #works on the 'panoramic_echo.jpg' and 'simple_echo.jpg' images in data folder
          #results panoramic_echo: 0.2631578947368421 0.2631578947368421
          #results simple_echo: 0.06097560975609756 0.06097560975609756
    """
    
    'Check if I is a 1-canal image'
    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    length, width = I.shape[0], I.shape[1]
    
    'Segmentation of the scale: cropping with empirical percentages and       '
    'binarization of the selection                                            '
    TCP = 0.1 #Top cropping percentage - #empirical percentage
    LCP = 0.5 #Left cropping percentage
    BCP = 0.7 #Bottom cropping percentage
    RCP = 0.1 #Right cropping percentage
    Scale_image = I[int(TCP * length):length-int(BCP * length),\
                    int(LCP * width):width-int(RCP * width)];                   
    Binar_I = cv2.threshold(Scale_image, 220., 255, cv2.THRESH_BINARY)[1]                
    
#    cv2.imshow('image',Binar_I);
#    cv2.waitKey(0) & 0xFF;
#    cv2.destroyAllWindows();    
    
    'Selection of the biggest dashes: contours of white objects are found as  '
    'well as minimal rectangles encapsulating each object. Conditions on the  '
    'size of these contours/bounding rectangles enable the removal of objects '
    'that are not the biggest dashes                                          '
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
   
    'Calculation of the minimal distances between two horizontal dashes and   '
    'two vertical dashes                                                      '
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

    'Factors to convert distance in pixels into distance in millimeters       '
    if horiz == 10000000.:
        calibFactorX = False
    else:
        calibFactorX = 10./horiz
    if vertic == 10000000.:
        calibFactorY = False
    else:
        calibFactorY = 10./vertic
        
    return calibFactorX, calibFactorY
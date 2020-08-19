from calibration.calib import autoCalibration
from preprocessing.cropping import autocropping
from preprocessing.preprocess import preprocessingApo
import aponeuroseslocation as apoL
import aponeuroses_contours as apoC
from MUFeaM import muscleThickness
import FaDe as FaDe

import cv2
import numpy as np
import tkinter.messagebox as tkbox

#################################EVENT FUNCTION################################

points = []
def _pickCoordinates(event):
    points.append((event.y,event.x))

########################PROCESSING OF SIMPLE US IMAGES#########################

#Open the image
#import aponeurosesdetection.data as apodat

#! RGBimage = apodat.simpleimg()
#RGBimage = cv2.imread('C:/Users/Antonio/Desktop/AponeurosesDetection/aponeurosesdetection/data/simple_echo.jpg', -1)

#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/01_Kevin/POST/architecture/Kevin_jamon_20180927_142814_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/31_romain/architecture/Romain_jamon_20181008_084409_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/fam_1/architecture/Julien_jamon_20180720_170728_image_bfs.jpg', -1)
RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/34_nicolas/architecture/Nicolas_post_20181210_105644_image_bfs.jpg', -1)
#RGBimage = data.simpleimg()
#################################################
cv2.imshow('Image to process', RGBimage)
process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
    Close all windows after clicking yes or no.', default = 'yes', icon='question')
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()

if process == True:
    #Calibrate the image
    calibX, calibY = autoCalibration(RGBimage)
    #################################################

    #Crop the image to keep essential US data
    USimage = autocropping(RGBimage, 10., 15., 12., 25., calibY, additionalCrop1 = 2., additionalCrop2 = 6.)

    #--------ask for validation; if cropping does not suit the user, ask the user for thresholds
    ok = False
    while ok == False:
        cv2.imshow('Cropped US image', USimage)
        ok = tkbox.askyesno('Need user validation', 'Do you validate the cropping ? If no, we will ask you new thresholds in the command prompt. After clicking yes or no, please close the image window to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()

        if ok == False:
            entry1 = input("Please enter an integer between 0 and 255 for minimum threshold to test for columns' mean. Recommended value is 10: ")
            entry2 = input("Please enter an integer between 0 and 255 for maximum threshold to test for columns' mean. Recommended value is 15: ")
            entry3 = input("Please enter an integer between 0 and 255 for minimum threshold to test for raws' mean. Recommended value is 12.: ")
            entry4 = input("Please enter an integer between 0 and 255 for maximum threshold to test for raws' mean. Recommended value is 25.: ")

            THRESH1 = int(entry1)
            THRESH2 = int(entry2)
            THRESH3 = int(entry3)
            THRESH4 = int(entry4)
            print(f'You entered the following thresholds: {THRESH1, THRESH2, THRESH3, THRESH4}')
            if THRESH1>255 or THRESH1<0 or THRESH2>255 or THRESH2<0 or THRESH3>255 or THRESH3<0 or\
                THRESH4>255 or THRESH4<0:
                raise ValueError('All thresholds must be integers between 0 and 255')
            USimage = autocropping(RGBimage, THRESH1, THRESH2, THRESH3, THRESH4, calibY)
    cv2.imwrite('C:/Users/Lisa Paillard/Desktop/nico644cropped.jpg', USimage)
    #################################################

    #Locate aponeuroses and find linear approximation of aponeuroses
    USimage_pp = preprocessingApo(USimage, 'localmean', 0, 41) #pre_processing
    aponeuroses_Linear, paramUp, paramLow, apoUp, apoLow = apoL.apoLocation(USimage_pp, thresh = None)
    
    USimage[apoUp[0], :,:] = [0,0,255]
    USimage[apoUp[1], :,:] = [0,0,255]
    USimage[apoLow[0], :,:] = [0,255,0]
    USimage[apoLow[1], :,:] = [0,255,0]
    cv2.imshow('interpolated aponeuroses', USimage)
    cv2.imshow('pp image', USimage_pp)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    #sub-images containing superficial and deep aponeurosis each
    upperApo = np.copy(USimage[apoUp[0]:apoUp[1],:])
    lowerApo = np.copy(USimage[apoLow[0]:apoLow[1],:])
    
    #Pre-processed sub-images of superficial and deep aponeuroses
    upperApo_pp = np.copy(USimage_pp[apoUp[0]:apoUp[1],:])
    lowerApo_pp = np.copy(USimage_pp[apoLow[0]:apoLow[1],:])
    
    #################################################
    
    # Get exact contour of each aponeurosis
    # First - we ask the user to click on points to get an initial contour
    #close to the real boundary of each aponeurosis (faster convergence)
    
    ini_upApo = apoC.initiateContour(upperApo_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1] - apoUp[0], 8])
    contourUp, nUp = apoC.activeContour(upperApo_pp,ini_upApo,0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
    print('Upper aponeurosis contour found in ', nUp, ' steps')
    if np.min(contourUp)<=0: #ask for validation of the contour, else linear approximation is used
        contourUpimage, contourPointsUp, obU = apoC.extractContour(contourUp, upperApo, offSetX = apoUp[0], offSetY = 0)
        cv2.imshow('Upper aponeurosis contour', contourUpimage)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contours ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        if valid == True: #B-spline to approximate aponeurosis if contour suits
            approxUp, coordUp = apoC.approximate(contourPointsUp, 'upper', upperApo_pp, d = 4)
        elif valid == False:
            yUp = np.arange(0,USimage.shape[1],1)
            xUp = np.int32(yUp*paramUp[0]+paramUp[1])
            coordUp = np.vstack((xUp, yUp)).T
    elif np.min(contourUp)>0: #use linear approximation because it means active contour model failed
        yUp = np.arange(0,USimage.shape[1],1)
        xUp = np.int32(yUp*paramUp[0]+paramUp[1])
        coordUp = np.vstack((xUp, yUp)).T
    
    ini_lowApo = apoC.initiateContour(lowerApo_pp, typeC = 'quadrangle_param', param = [paramLow[0], paramLow[1] - apoLow[0], 8])
    contourLow, nLow = \
        apoC.activeContour(lowerApo_pp, ini_lowApo, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
    print('Lower aponeurosis contour found in ', nLow, ' steps')
    if np.min(contourLow)<=0: #ask for validation of the contour, else linear approximation is used
        contourLowimage, contourPointsLow, obL = apoC.extractContour(contourLow, lowerApo, offSetX = apoLow[0], offSetY = 0)
        cv2.imshow('Lower aponeurosis contour', contourLowimage)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contours ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()    
        if valid == True: #B-spline to approximate aponeurosis if contour suits
            approxLow, coordLow = apoC.approximate(contourPointsLow, 'lower', lowerApo_pp, d = 4)
        elif valid == False:
            yLow = np.arange(0,USimage.shape[1],1)
            xLow = np.int32(yLow*paramLow[0]+paramLow[1])
            coordLow = np.vstack((xLow, yLow)).T
    elif np.min(contourLow)>0: #use linear approximation because it means active contour model failed
        yLow = np.arange(0,USimage.shape[1],1)
        xLow = np.int32(yLow*paramLow[0]+paramLow[1])
        coordLow = np.vstack((xLow, yLow)).T        
    
    #Muscle thickness calculation
    a, thickness, spline_thickness = muscleThickness(USimage, coordUp, coordLow, 0, USimage.shape[1]-1, calibX, calibY)
    
    #visualization:
    for index in range(coordUp.shape[0]):
        USimage[coordUp[index][0], coordUp[index][1],:] = [0,0,255]
        USimage[coordLow[index][0], coordLow[index][1], :] = [0,0,255]
    cv2.imshow('interpolated aponeuroses', USimage)
    import matplotlib.pyplot as plt
    plt.plot(a, thickness)
    plt.show()
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    #Crop to focus in between aponeuroses
    #Erase information above superficial apo and below deep apo
    
    for col in range(USimage.shape[1]):
        line = 0
        while line != coordUp[col][0]:
            USimage[line, col, :] = [0,0,0]
            line = line + 1
        line = USimage.shape[0] - 1
        while line != coordLow[col][0]:
            USimage[line, col, :] = [0,0,0]
            line = line - 1
    crop1 = np.amax(coordUp[:,0]) + 1 
    crop2 = np.amin(coordLow[:,0])
    USimage = USimage[crop1:crop2, :]
    cv2.imshow('Region of interest', USimage)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    #Enhance tube-like structures with MVEF method - Frangi
    #considering that fascicle diameter is between 0.15 mm and 0.5 mm,
    #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
    sca = np.arange(round(0.3/calibX*2)/2, round(0.5/calibX*2)/2, 0.5)
    print('MVEF running')
    #
    imgMVEF = FaDe.MVEF_2D(255-USimage, sca, [0.5, 0])
    #
    #threshold
    threshMVEF_percent = 85
    threshMVEF = np.percentile(imgMVEF, threshMVEF_percent)
    I2 = cv2.threshold(imgMVEF, threshMVEF, 255, cv2.THRESH_BINARY)[1]
    
    #locate muscle fascicles
    fasc, Nsnippets = FaDe.locateFascicle(I2[:I2.shape[0],:], calibX, calibY, USimage)
    
    for i in range(len(fasc)):
        if len(Nsnippets[i])>1:
            for j in range(fasc[i].shape[0]):
                USimage[fasc[i][j,1],fasc[i][j,0],:] = [255,255,255]
        else:
            for j in range(fasc[i].shape[0]):
                USimage[fasc[i][j,1],fasc[i][j,0],:] = [0,255,0]
    
    cv2.imshow('final', I2)
    cv2.imshow('img', USimage)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    #cv2.imshow('MVEF', imgMVEF)
    #cv2.imshow('final', I2)
    #cv2.imshow('img', USimage)
    #cv2.waitKey(0) & 0xFF
    #cv2.destroyAllWindows()
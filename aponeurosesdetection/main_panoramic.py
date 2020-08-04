from calibration.calib import autoCalibration
from preprocessing.cropping import autocropping, manualcropping
from preprocessing.preprocessing import preprocessingApo
import aponeuroseslocation as apoL
import aponeuroses_contours as apoC

import cv2
import numpy as np
import tkinter.messagebox as tkbox
import tkinter as tk
from PIL import ImageTk, Image

#################################EVENT FUNCTION################################

points = []
def _pickCoordinates(event):
    points.append((event.y,event.x))

#######################PROCESSING OF PANORAMIC US IMAGES#######################

#Open the image
RGBimage_p = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/code/semi manuel/Romain_jamon_20181008_084851_image_bfp.jpg', -1)

#################################################

#Calibrate the image
calibX_p, calibY_p = autoCalibration(RGBimage_p)

#################################################

#Crop the image thanks to manual preprocessing
#(textfile containing points coordinates)
USimage_p = manualcropping(RGBimage_p, 'C:/Users/Lisa Paillard/Desktop/Pour TFE/code/semi manuel/Romain_jamon_20181008_084851_image_bfp.txt')

#################################################

#resize
initialsize = (USimage_p.shape[1], USimage_p.shape[0])
PERCENTAGE = 160
newWidth = int(USimage_p.shape[1]*PERCENTAGE/100)
newHeight = int(USimage_p.shape[0]*PERCENTAGE/100)
USimage_p = cv2.resize(src = USimage_p, dsize = (newWidth, newHeight), interpolation = cv2.INTER_CUBIC)
calibX_p = calibX_p * PERCENTAGE / 100
calibY_p = calibY_p * PERCENTAGE / 100

#################################################

#Sampling to analyze band-by-band the image on MAXBANDS/NBANDS of the width
NBANDS = 6 #number of bands used to sample the image
MAXBAND = 3 #last band to be analyzed (start counting from the left)
            #all bands cannot be analyzed because deep aponeurosis is not visible enough

sampleSize = int(USimage_p.shape[1] / NBANDS)
        # sample1 = USimage_p[:,:sampleSize]
        # sample2 = USimage_p[:,sampleSize:2*sampleSize]
        # etc

#################################################

#Preprocessing, location of aponeuroses and linear approximation of aponeuroses

splines = np.zeros((MAXBAND, 2, 2, sampleSize))
for i in range(MAXBAND):
    #
    apo_lin, paramUp, paramDeep, locUp, locDeep = apoL.apoLocation(USimage_p[:, i*sampleSize:(i+1)*sampleSize], 200.)

    #
    UAi = np.copy(USimage_p[locUp[0]:locUp[1], i*sampleSize:(i+1)*sampleSize]) #upper aponeurosis in sample i
    UAi_pp = preprocessingApo(UAi, 'localmidgrey', 0, 41)
    points = []
    window = tk.Tk()
    canvas = tk.Canvas(window, width = UAi_pp.shape[1], height = UAi_pp.shape[0])      
    canvas.pack()      
    displayI = ImageTk.PhotoImage(image = Image.fromarray(UAi_pp), master = window)     
    canvas.create_image(0,0, anchor='nw', image=displayI) 
    window.bind('<Button-1>',_pickCoordinates) 
    window.mainloop()
    points = np.array(points)
    print(points)
    iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'set_of_points', setPoints = points)
    contourUp_i, nUp_i = \
        apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
    contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi)
    print('Upper aponeurosis contour found in ', nUp_i, ' steps')

    #
    DAi = np.copy(USimage_p[locDeep[0]:locDeep[1], i*sampleSize:(i+1)*sampleSize]) #deep aponeurosis in sample i
    DAi_pp = preprocessingApo(DAi, 'localmidgrey', 0, 41)
    points = []
    window = tk.Tk()
    canvas = tk.Canvas(window, width = DAi_pp.shape[1], height = DAi_pp.shape[0])      
    canvas.pack()      
    displayI = ImageTk.PhotoImage(image = Image.fromarray(DAi_pp), master = window)     
    canvas.create_image(0,0, anchor='nw', image=displayI) 
    window.bind('<Button-1>',_pickCoordinates) 
    window.mainloop()
    points = np.array(points)
    print(points)
    iniDeep_i = apoC.initiateContour(DAi_pp, typeC = 'set_of_points', setPoints = points)
    contourDeep_i, nDeep_i = \
        apoC.activeContour(DAi_pp, iniDeep_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
    contourDeepimage_i, contourPointsDeep_i = apoC.extractContour(contourDeep_i, DAi)
    print('Deep aponeurosis contour found in ', nDeep_i, ' steps')

    #validation of the contours
    #if not validated, linear approximation is used
    cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
    cv2.imshow('Lower aponeurosis contour on sample i', contourDeepimage_i)
    valid = tkbox.askyesno('Need user validation', 'Do you validate the contours ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()

    if valid == True:
        approxUp_i, xUp_i, yUp_i = apoC.approximate(contourPointsUp_i, 'upper', UAi_pp, d = 4)
        approxDeep_i, xDeep_i, yDeep_i = apoC.approximate(contourPointsDeep_i, 'lower', DAi_pp, d = 4)
        #transform back to USimage coordinates:
        xUp_i = xUp_i + locUp[0]
        xDeep_i = xDeep_i + locDeep[0]

    elif valid == False:
        yUp_i = np.arange(0,UAi.shape[1]-1,1)
        xUp_i = np.int32(yUp_i*paramUp[0]+paramUp[1])
        yDeep_i = np.arange(0,DAi.shape[1]-1,1)
        xDeep_i = np.int32(yDeep_i*paramDeep[0]+paramDeep[1])

    splines[i, 0, 0, :] = xUp_i
    splines[i, 0, 1, :] = yUp_i
    splines[i, 1, 0, :] = xDeep_i
    splines[i, 1, 1, :] = yDeep_i

    #Visualization
    for index in range(yUp_i.shape[0]):
        USimage_p[xUp_i[index], yUp_i[index] + i * sampleSize, :] = [0, 255, 0]
        USimage_p[xDeep_i[index], yDeep_i[index] + i * sampleSize, :] = [0, 255, 0]

#NOT FINISHED : MISSING LAST PART TO INTERPOLATE ENTIRE APONEUROSES 


cv2.imshow('interpolated aponeuroses', USimage_p)
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()
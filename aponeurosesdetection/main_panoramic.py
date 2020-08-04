from calibration.calib import autoCalibration
from preprocessing.cropping import autocropping, manualcropping
from preprocessing.preprocessing import preprocessingApo
import aponeurosesdetection.aponeuroseslocation as apoL
import aponeurosesdetection.aponeuroses_contours as apoC

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

#Sampling to analyze band-by-band the image on 2/3 of the width
sampleSize = int(USimage_p.shape[1] / 6)
# sample1 = USimage_p[:,:sampleSize]
# sample2 = USimage_p[:,sampleSize:2*sampleSize]
# sample3 = USimage_p[:,2*sampleSize:3*sampleSize]
# sample4 = USimage_p[:,3*sampleSize:4*sampleSize]
# sample5 = USimage_p[:,4*sampleSize:5*sampleSize]
# sample6 = USimage_p[:,5*sampleSize:USimage_p.shape[1]]

#################################################

#Preprocessing, location of aponeuroses and linear approximation of aponeuroses
locations = np.zeros((4, 2))
linesParam = np.zeros((4, 2))
for i in range(4):
    #sampling and preprocessing
    ppSample_i = preprocessingApo(USimage_p[:, i*sampleSize:(i+1)*sampleSize], 'localmidgrey', 0, 41)
    apo_lin, paramUp, paramDeep, locUp, locDeep = apoL.apoLocation(ppSample_i, 200.)
    locations[i, 0], locations[i, 1] = locUp, locDeep
    linesParam[i, 0], linesParam[i, 1] = paramUp, paramDeep

    points = []
    window = tk.Tk()
    canvas = tk.Canvas(window, width = ppSample_i.shape[1], height = ppSample_i.shape[0])      
    canvas.pack()      
    displayI = ImageTk.PhotoImage(image = Image.fromarray(ppSample_i), master = window)     
    canvas.create_image(0,0, anchor='nw', image=displayI) 
    window.bind('<Button-1>',_pickCoordinates) 
    window.mainloop()
    points = np.array(points)
    print(points)
    iniUpper_i = apoC.initiateContour(ppSample_i, typeC = 'set_of_points', setPoints = points)
    contourUp_i, nUp_i = \
        apoC.activeContour(ppSample_i, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
    contourUpimage_i, contourPoints_i = apoC.extractContour(contourUp_i, USimage_p[:, i*sampleSize:(i+1)*sampleSize])
    print('Lower aponeurosis contour found in ', nUp_p, ' steps')

"""
#################################################

window = tk.Tk()
canvas = tk.Canvas(window, width = sample1.shape[1], height = sample1.shape[0])      
canvas.pack()      
displayI = ImageTk.PhotoImage(image = Image.fromarray(sample1p), master = window)     
canvas.create_image(0,0, anchor='nw', image=displayI) 
window.bind('<Button-1>',_pickCoordinates) 
window.mainloop()
points = np.array(points)
print(points)
iniUpper_p = apoC.initiateContour(sample1p, typeC = 'set_of_points', setPoints = points)
contourUp_p, nUp_p = \
    apoC.activeContour(sample1p, iniUpper_p, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
contourUpimage_p, contourPoints_p = apoC.extractContour(contourUp_p, sample1)
print('Lower aponeurosis contour found in ', nUp_p, ' steps')

cv2.imshow('contour up image', contourUpimage_p)

cv2.imshow('Manually cropped panoramic image', USimage_p)
#cv2.imshow('Pre-processed image', imagePrepro_p)

cv2.imshow('preprocessed', sample1p)
cv2.imshow('contour', contourUpimage_p)
# cv2.imshow('sample2', sample2)
# cv2.imshow('sample3', sample3)
# cv2.imshow('sample4', sample4)
# cv2.imshow('sample5', sample5)
# cv2.imshow('sample6', sample6)
"""
cv2.imshow('sample1', sample1)
cv2.imshow('apoup',upperApo)
cv2.imshow('deep apo',lowerApo)


cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()
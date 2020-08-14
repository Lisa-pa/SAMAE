from aponeurosesdetection.calibration.calib import autoCalibration
from aponeurosesdetection.preprocessing.cropping import autocropping
from aponeurosesdetection.preprocessing.preprocess import preprocessingApo
import aponeurosesdetection.aponeuroseslocation as apoL
import aponeurosesdetection.aponeuroses_contours as apoC
from aponeurosesdetection.MUFeaM import muscleThickness
import aponeurosesdetection.FaDe as FaDe

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

RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/AponeurosesDetection/aponeurosesdetection/data/simple_echo.jpg', -1)
#RGBimage = data.simpleimg()
#################################################

#Calibrate the image
calibX, calibY = autoCalibration(RGBimage)

#################################################

#Crop the image to keep essential US data
USimage = autocropping(RGBimage, 10., 15., 12., 25., calibY)

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

#################################################

#Locate aponeuroses and find linear approximation of aponeuroses
aponeuroses_Linear, paramUp, paramLow, apoUp, apoLow = apoL.apoLocation(USimage, 220.)
upperApo = np.copy(USimage[apoUp[0]:apoUp[1],:])
lowerApo = np.copy(USimage[apoLow[0]:apoLow[1],:])

#################################################

#Pre-process the sub-images of upper and lower aponeuroses
upperApo_pp = preprocessingApo(upperApo, 'localmean', 0, 41)
lowerApo_pp = preprocessingApo(lowerApo, 'localmean', 0, 41)

#################################################

# Get exact contour of each aponeurosis
# First - we ask the user to click on points to get an initial contour
#close to the real boundary of each aponeurosis (faster convergence)

ini_upApo = apoC.initiateContour(upperApo_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1] - apoUp[0], 8])
contourUp, nUp = apoC.activeContour(upperApo_pp,ini_upApo,0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
contourUpimage, contourPointsUp, obU = apoC.extractContour(contourUp, upperApo, offSetX = apoUp[0], offSetY = 0)
print('Upper aponeurosis contour found in ', nUp, ' steps')

ini_lowApo = apoC.initiateContour(lowerApo_pp, typeC = 'quadrangle_param', param = [paramLow[0], paramLow[1] - apoLow[0], 8])
contourLow, nLow = \
    apoC.activeContour(lowerApo_pp, ini_lowApo, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
contourLowimage, contourPointsLow, obL = apoC.extractContour(contourLow, lowerApo, offSetX = apoLow[0], offSetY = 0)
print('Lower aponeurosis contour found in ', nLow, ' steps')

#################################################

#validation of the contours
#if not validated, linear approximation is used
cv2.imshow('Upper aponeurosis contour', contourUpimage)
cv2.imshow('Lower aponeurosis contour', contourLowimage)
valid = tkbox.askyesno('Need user validation', 'Do you validate the contours ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()

#################################################

#B-spline to approximate each aponeurosis if contours suit

if valid == True:
    approxUp, coordUp = apoC.approximate(contourPointsUp, 'upper', upperApo_pp, d = 4)
    approxLow, coordLow = apoC.approximate(contourPointsLow, 'lower', lowerApo_pp, d = 4)

elif valid == False:
    yUp = np.arange(0,USimage.shape[1]-1,1)
    xUp = np.int32(yUp*paramUp[0]+paramUp[1])
    coordUp = np.vstack((xUp, yUp)).T
    yLow = np.arange(0,USimage.shape[1]-1,1)
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
crop1 = np.amin(coordUp[:,0])
crop2 = np.amax(coordLow[:,0])
USimage = USimage[crop1:crop2, :]
cv2.imshow('Erased', USimage)
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()

#Enhance tube-like structures with MVEF method - Frangi
#considering that fascicle diameter is between 0.15 mm and 0.5 mm,
#the following list is the equivalent interval in pixels, with a step of 0.5 pixel
sca = np.arange(round(0.3/calibX*2)/2, round(0.5/calibX*2)/2, 0.5)
print(sca)
#
imgMVEF = FaDe.MVEF_2D(255-USimage, sca, [0.5, 0])
cv2.imshow('MVEF', imgMVEF)
cv2.imshow('img', USimage)
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()
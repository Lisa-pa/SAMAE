from aponeurosesdetection.calibration.calib import autoCalibration
from aponeurosesdetection.preprocessing.cropping import manualcropping
from aponeurosesdetection.preprocessing.preprocess import preprocessingApo
import aponeurosesdetection.aponeuroseslocation as apoL
import aponeurosesdetection.aponeuroses_contours as apoC
from aponeurosesdetection.MUFeaM import muscleThickness

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

#Crop the image thanks to manual labelling
#(need text file containing points coordinates)
#pt_intersection is the point where aponeuroses meet
USimage_p, pt_intersection = manualcropping(RGBimage_p, 'C:/Users/Lisa Paillard/Desktop/Pour TFE/code/semi manuel/Romain_jamon_20181008_084851_image_bfp.txt')

#################################################

#resize and update calibration factors and pt_intersection's coordinates
initialsize = (USimage_p.shape[1], USimage_p.shape[0])
PERCENTAGE = 160
newWidth = int(initialsize[0]*PERCENTAGE/100)
newHeight = int(initialsize[1]*PERCENTAGE/100)
USimage_p = cv2.resize(src = USimage_p, dsize = (newWidth, newHeight), interpolation = cv2.INTER_CUBIC)

calibX_p = calibX_p / PERCENTAGE * 100
calibY_p = calibY_p / PERCENTAGE * 100

pt_intersection = (pt_intersection[0] * PERCENTAGE / 100, pt_intersection[1] * PERCENTAGE / 100)

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

contoursUp = []
contoursDeep = []

for i in range(NBANDS):
    #
    #
    if i < MAXBAND:
        #find where each aponeurosis is
        apo_lin, paramUp, paramDeep, locUp, locDeep = apoL.apoLocation(USimage_p[:, i*sampleSize:(i+1)*sampleSize], 200.)
        
        #sub-images of each aponeurosis
        UAi = np.copy(USimage_p[locUp[0]:locUp[1], i*sampleSize:(i+1)*sampleSize]) #upper aponeurosis in sample i
        UAi_pp = preprocessingApo(UAi, 'localmidgrey', 0, 41)
        DAi = np.copy(USimage_p[locDeep[0]:locDeep[1], i*sampleSize:(i+1)*sampleSize]) #deep aponeurosis in sample i
        DAi_pp = preprocessingApo(DAi, 'localmidgrey', 0, 41)

        ### upper apo
        #ask for initial contour
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
        #Calculate contour with active contour model
        iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'set_of_points', setPoints = points)
        contourUp_i, nUp_i = \
            apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
        contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = locUp[0], offSetY = i * sampleSize)
        print('Upper aponeurosis contour found in ', nUp_i, ' steps')
        #ask for manual validation of the contour
        cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        if valid == True:
            #add contour_i to list 'contoursUp'
            for elem in contourPointsUp_i:
                contoursUp.append(elem)

        ### deep apo
        #ask for initial contour
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
        #calculate contour with active contour model
        iniDeep_i = apoC.initiateContour(DAi_pp, typeC = 'set_of_points', setPoints = points)
        contourDeep_i, nDeep_i = \
            apoC.activeContour(DAi_pp, iniDeep_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
        contourDeepimage_i, contourPointsDeep_i = apoC.extractContour(contourDeep_i, DAi, offSetX = locDeep[0], offSetY = i * sampleSize)
        print('Deep aponeurosis contour found in ', nDeep_i, ' steps')
        #ask for manual validation of the contour
        cv2.imshow('Lower aponeurosis contour on sample i', contourDeepimage_i)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contours ? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        if valid == True:
            #add contour_i to contourDeep
            for elem in contourPointsDeep_i:
                contoursDeep.append(elem)
    else:
        ### we only consider upper aponeurosis, and we know that it stands in the first half of the image
        UAi = np.copy(USimage_p[:int(USimage_p.shape[0]/2), i*sampleSize:min((i+1)*sampleSize, int(pt_intersection[1]))]) #upper aponeurosis in sample i
        UAi_pp = preprocessingApo(UAi, 'localmidgrey', 0, 41)
        #ask for initial contour
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
        #calculate contour with active contour model
        iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'set_of_points', setPoints = points)
        contourUp_i, nUp_i = \
            apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
        contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = 0, offSetY = i * sampleSize)
        print('Upper aponeurosis contour found in ', nUp_i, ' steps')
        #ask for manual validation of the contour
        cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        if valid == True:
            #add contour_i to list 'contoursUp'
            for elem in contourPointsUp_i:
                contoursUp.append(elem)

contoursUp.append(pt_intersection)
contoursDeep.append(pt_intersection)

spline_up, coordUp = apoC.approximate(contoursUp, 'upper', USimage_p, d = 3)
spline_deep, coordDeep = apoC.approximate(contoursDeep, 'lower', USimage_p, d = 3)

#muscle thickness measurement
a, thickness, thickness_spline = muscleThickness(USimage_p, coordUp, coordDeep, 0, int(pt_intersection[1]), calibX_p, calibY_p)

#Visualization
for index in range(coordUp.shape[0]):
    USimage_p[coordUp[index][0], coordUp[index][1], :] = [0, 255, 0]
    USimage_p[coordDeep[index][0], coordDeep[index][1], :] = [0, 255, 0]
cv2.imshow('interpolated aponeuroses', USimage_p)
import matplotlib.pyplot as plt
plt.plot(a, thickness)
plt.show()

cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()
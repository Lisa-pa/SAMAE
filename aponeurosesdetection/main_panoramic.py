from calibration.calib import autoCalibration
from preprocessing.cropping import manualcropping
from preprocessing.preprocess import preprocessingApo
import aponeuroseslocation as apoL
import aponeuroses_contours as apoC
import MUFeaM
import FaDe

import cv2
import numpy as np
import tkinter.messagebox as tkbox
#import tkinter as tk
#from PIL import ImageTk, Image

#################################EVENT FUNCTION################################

points = []
def _pickCoordinates(event):
    points.append((event.y,event.x))

#######################PROCESSING OF PANORAMIC US IMAGES#######################

#Open the image
RGBimageP = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/code/semi manuel/Romain_jamon_20181008_084851_image_bfp.jpg', -1)

cv2.imshow('Image to process', RGBimageP)
process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
    Close all windows after clicking yes or no.', default = 'yes', icon='question')
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()

if process == True:
    
    #################################################
    
    #Calibrate the image
    calibX, calibY = autoCalibration(RGBimageP)
    
    #################################################
    
    #Crop the image thanks to manual labelling
    #(need text file containing points coordinates)
    #pt_intersection is the point where aponeuroses meet
    USimageP, pt_intersection = manualcropping(RGBimageP, 'C:/Users/Lisa Paillard/Desktop/Pour TFE/code/semi manuel/Romain_jamon_20181008_084851_image_bfp.txt')
    
    #################################################
    
    #resize and update calibration factors and pt_intersection's coordinates
    initialsize = (USimageP.shape[1], USimageP.shape[0])
    PERCENTAGE = 160
    newWidth = int(initialsize[0]*PERCENTAGE/100)
    newHeight = int(initialsize[1]*PERCENTAGE/100)
    USimageP = cv2.resize(src = USimageP, dsize = (newWidth, newHeight), interpolation = cv2.INTER_CUBIC)
    
    calibX = calibX / PERCENTAGE * 100
    calibY = calibY / PERCENTAGE * 100
    
    pt_intersection = (pt_intersection[0] * PERCENTAGE / 100, pt_intersection[1] * PERCENTAGE / 100)
    
    #################################################
    
    #Sampling to analyze band-by-band the image on MAXBANDS/NBANDS of the width
    NBANDS = 6 #number of bands used to sample the image
    MAXBAND = 3 #last band to be analyzed (start counting from the left)
                #all bands cannot be analyzed because deep aponeurosis is not visible enough
    sampleSize = int(USimageP.shape[1] / NBANDS)
            # sample1 = USimageP[:,:sampleSize]
            # sample2 = USimageP[:,sampleSize:2*sampleSize]
            # etc
    
    #################################################
    
    #Preprocessing
    USimageP_pp = preprocessingApo(USimageP, 'localmidgrey', 0, 41)
    
    #################################################   
    #location of aponeuroses and linear approximation of aponeuroses
    
    contoursUp = []
    contoursDeep = []
    
    for i in range(NBANDS):
        #
        #
        if i < MAXBAND:
            #find where each aponeurosis is
            paramUp, paramDeep, locUp, locDeep = apoL.twoApoLocation(USimageP_pp[:, i*sampleSize:(i+1)*sampleSize], None)
            
            ### upper apo
            #sub-images of superficial/Up aponeurosis
            UAi = np.copy(USimageP[locUp[0]:locUp[1], i*sampleSize:(i+1)*sampleSize])
            UAi_pp = np.copy(USimageP_pp[locUp[0]:locUp[1], i*sampleSize:(i+1)*sampleSize]) #upper aponeurosis in sample i
    
            #Calculate contour with active contour model
            iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1]-locUp[0], 10])       
            
            contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
            print('Upper aponeurosis contour found in ', nUp_i, ' steps')


            if np.amin(contourUp_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1] - locUp[0], 40])
                contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                print('Upper aponeurosis contour found in ', nUp_i, ' steps')

            if np.amin(contourUp_i)<=0: #if the contour exists, extract it
                contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = locUp[0], offSetY = i * sampleSize)

                #ask for manual validation of the contour
                cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
                valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                cv2.waitKey(0) & 0xFF
                cv2.destroyAllWindows()
                
                if valid == False:
                    #try a second time with a new initial contour
                    points = np.array([[0, 0], [int(UAi_pp.shape[0]/2), 0], [int(UAi_pp.shape[0]/2), UAi_pp.shape[1]], [0, UAi_pp.shape[1]]])
                    iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'set_of_points', setPoints = points)
                    contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Upper aponeurosis contour found in ', nUp_i, ' steps')
                    if np.amin(contourUp_i) <= 0 :
                        contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = locUp[0], offSetY = i * sampleSize)
                        #ask for manual validation of the contour
                        cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
                        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                        cv2.waitKey(0) & 0xFF
                        cv2.destroyAllWindows()                
                
                
                if valid == True:
                    #add contour_i to list 'contoursUp'
                    for elem in contourPointsUp_i:
                        contoursUp.append(elem)
    
            ### deep apo
            #sub-images of deep aponeurosis
            DAi = np.copy(USimageP[locDeep[0]:locDeep[1], i*sampleSize:(i+1)*sampleSize])
            DAi_pp = np.copy(USimageP_pp[locDeep[0]:locDeep[1], i*sampleSize:(i+1)*sampleSize]) #deep aponeurosis in sample i
            
            #calculate contour with active contour model
            iniDeep_i = apoC.initiateContour(DAi_pp, typeC = 'quadrangle_param', param = [paramDeep[0], paramDeep[1]-locDeep[0], 10])

            contourDeep_i, nDeep_i = apoC.activeContour(DAi_pp, iniDeep_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
            print('Deep aponeurosis contour found in ', nDeep_i, ' steps')
            
            if np.amin(contourDeep_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                iniDeep_i = apoC.initiateContour(DAi_pp, typeC = 'quadrangle_param', param = [paramDeep[0], paramDeep[1] - locDeep[0], 40])
                contourDeep_i, nDeep_i = apoC.activeContour(DAi_pp, iniDeep_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                print('Upper aponeurosis contour found in ', nDeep_i, ' steps')

            if np.amin(contourDeep_i)<=0: #if the contour exists, extract it
                contourDeepimage_i, contourPointsDeep_i = apoC.extractContour(contourDeep_i, DAi, offSetX = locDeep[0], offSetY = i * sampleSize)
            
                #ask for manual validation of the contour
                cv2.imshow('Lower aponeurosis contour on sample i', contourDeepimage_i)
                valid = tkbox.askyesno('Need user validation', 'Do you validate the contours? After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                cv2.waitKey(0) & 0xFF
                cv2.destroyAllWindows()
                
                if valid == False:
                    #try a second time with a new initial contour
                    points = np.array([[int(DAi_pp.shape[0]/2), 0], [DAi_pp.shape[0], 0], [DAi_pp.shape[0],DAi_pp.shape[1]], [int(DAi_pp.shape[0]/2), DAi_pp.shape[1]]])
                    iniDeep_i = apoC.initiateContour(DAi_pp, typeC = 'set_of_points', setPoints = points)
                    contourDeep_i, nDeep_i = apoC.activeContour(DAi_pp, iniDeep_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Deep aponeurosis contour found in ', nDeep_i, ' steps')
                    if np.amin(contourDeep_i) <= 0 :
                        contourDeepimage_i, contourPointsDeep_i = apoC.extractContour(contourDeep_i, DAi, offSetX = locDeep[0], offSetY = i * sampleSize)
                        #ask for manual validation of the contour
                        cv2.imshow('Deep aponeurosis contour on sample i', contourDeepimage_i)
                        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                        cv2.waitKey(0) & 0xFF
                        cv2.destroyAllWindows()                
                
                if valid == True:
                    #add contour_i to contourDeep
                    for elem in contourPointsDeep_i:
                        contoursDeep.append(elem)
            
        else:
            ### we only consider upper aponeurosis, and we know that it stands in the first half of the image
            UAi = np.copy(USimageP[:int(USimageP.shape[0]/2), i*sampleSize:min((i+1)*sampleSize, int(pt_intersection[1]))])
            UAi_pp = np.copy(USimageP_pp[:int(USimageP.shape[0]/2), i*sampleSize:min((i+1)*sampleSize, int(pt_intersection[1]))])
            
            #initiate contour
            param, loc = apoL.oneApoLocation(UAi_pp)
            iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'quadrangle_param', param = [param[0], param[1], 10])
                        
            #calculate contour with active contour model
            contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
            print('Upper aponeurosis contour found in ', nUp_i, ' steps')
            
            if np.amin(contourUp_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'quadrangle_param', param = [param[0], param[1], 40])
                contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                print('Upper aponeurosis contour found in ', nUp_i, ' steps')
                  
            if np.amin(contourUp_i) <= 0 :
                contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = 0, offSetY = i * sampleSize)
                
                #ask for manual validation of the contour
                cv2.imshow('Upper aponeurosis contour on sample i', contourUpimage_i)
                valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                cv2.waitKey(0) & 0xFF
                cv2.destroyAllWindows()
                
                if valid == False:
                    #try a second time with a new initial contour
                    points = np.array([[0, 0], [int(UAi_pp.shape[0]/2), 0],[int(UAi_pp.shape[0]/2), UAi_pp.shape[1]], [0, UAi_pp.shape[1]]])
                    iniUpper_i = apoC.initiateContour(UAi_pp, typeC = 'set_of_points', setPoints = points)
                    contourUp_i, nUp_i = apoC.activeContour(UAi_pp, iniUpper_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Upper aponeurosis contour found in ', nUp_i, ' steps')
                    if np.amin(contourUp_i) <= 0 :
                        contourUpimage_i, contourPointsUp_i = apoC.extractContour(contourUp_i, UAi, offSetX = 0, offSetY = i * sampleSize)
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
    
    #Interpolation and extrapolation
    spline_up = apoC.approximateApo(contoursUp, 'upper', USimageP, d = 1)
    spline_deep = apoC.approximateApo(contoursDeep, 'lower', USimageP, d = 1)
    
    
    #muscle thickness measurement
    abscissa, thickness, thickness_spline = MUFeaM.muscleThickness(I = USimageP, spl1 = spline_up, spl2 = spline_deep, start = 0, end = int(pt_intersection[1]), calibV = calibX, calibH = calibY)


    #Calculate coordinates of aponeuroses
    coordUp = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_up, interval = [0, int(pt_intersection[1])])
    coordDeep = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_deep, interval = [0, int(pt_intersection[1])])
    
    
   #Fascicles detection
    all_snippets = []
    all_snippets_line = []
    
    for i in range(NBANDS-1):
        
        minRow = np.amax(coordUp[i*sampleSize:(i+1)*sampleSize,0])
        maxRow = np.amin(coordDeep[i*sampleSize:(i+1)*sampleSize,0])
        ROI = USimageP[minRow:maxRow, i*sampleSize:min((i+1)*sampleSize, pt_intersection[1]), :]
        
        #Enhance tube-like structures with MVEF method - Frangi -
        #let's consider that fascicle diameter is between 0.3 mm and 0.5 mm,
        #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
        print('MVEF running')
        sca = np.arange(round(0.3/calibX), round(0.5/calibX), 0.5)
        MVEF_image = FaDe.MVEF_2D(255-ROI, sca, [0.5, 0])
        #
        #threshold
        threshMVEF_percent = 85
        threshMVEF = np.percentile(MVEF_image, threshMVEF_percent)
        MVEF_image2 = cv2.threshold(MVEF_image, threshMVEF, 255, cv2.THRESH_BINARY)[1]
      
        #locate muscle snippets and filter them
        snippets, snippets_line = FaDe.locateSnippets(MVEF_image2, calibX, calibY,\
                                                      minLength = 4, rangeAngles = [6,40],\
                                                      percentageAlign = 0.80,\
                                                      offSetX = minRow, offSetY = i*sampleSize)
        all_snippets = all_snippets + snippets
        all_snippets_line = all_snippets_line + snippets_line

    fasc = FaDe.combineSnippets(USimageP, all_snippets, all_snippets_line, thresh_alignment = 5, thresh_length = 1000)

      
    #transform the snippets, which are in fact the contours of the fascicles,
    #into lines by taking the mean of the contour, that is the mean on each
    #column of the contour. When branches exist, the mean is not computed (the
    #region is ignored)

    averages = [] 
    for i in range(len(fasc)):
        averages.append(FaDe.contourAverage(fasc[i]))
        
    #interpolations to get fascicles' curve
    splines_fasc = FaDe.approximateFasc(I = USimageP, typeapprox = 'Bspline', listF = averages, d = 1)

    #Location of fascicles: (in mm from aponeuroses intersection point)
    loc_fasc = MUFeaM.locateFasc(splines_fasc, pt_intersection, calibY)
    print(loc_fasc)

    #intersections of fascicles with aponeuroses (in pixels)
    intersecL = MUFeaM.findIntersections(spline_deep, splines_fasc, USimageP, typeI = 'panoramic')
    intersecU = MUFeaM.findIntersections(spline_up, splines_fasc, USimageP, typeI = 'panoramic')
    
    #Pennation angles calculation (in degree)
    PA_up = MUFeaM.pennationAngles(spline_up, splines_fasc, intersecU, calibX, calibY)
    PA_low = MUFeaM.pennationAngles(spline_deep, splines_fasc, intersecL, calibX, calibY)
    
    #Fascicles length calculation (in mm)
    fasc_length = MUFeaM.fasciclesLength(splines_fasc, intersecU, intersecL, calibX, calibY)

    #Tuples containing info per fascicle
    muscle_fascicles = []
    for index in range(len(splines_fasc)):
        typeIU = 'out of image'
        typeIL = 'out of image'
        typeFL = 'out of image'
        if intersecU[index][0] >= 0 and\
            intersecU[index][0] < USimageP.shape[0] and\
            intersecU[index][1] >= 0 and\
            intersecU[index][1] < USimageP.shape[1]:
            typeIU = 'in image'
        if intersecL[index][0] >= 0 and\
            intersecL[index][0] < USimageP.shape[0] and\
            intersecL[index][1] >= 0 and\
            intersecL[index][1] < USimageP.shape[1]:
            typeIL = 'in image'
        if typeIU == 'in image' and typeIL == 'in image':
            typeFL = 'in image'

        tup = (index, loc_fasc[index], typeIU, intersecU[index], PA_up[index], typeIL, intersecL[index], PA_low[index], typeFL, fasc_length[index]) 
        muscle_fascicles.append(tup)


















    #####
    #Visualization
    for index in range(coordUp.shape[0]):
        if coordUp[index][0] >= 0 and coordUp[index][0] < USimageP.shape[0]:
            USimageP[coordUp[index][0], coordUp[index][1], :] = [0, 255, 0]
            if coordUp[index][0] +1 >= 0 and coordUp[index][0]+1 < USimageP.shape[0]:
                USimageP[coordUp[index][0]+1, coordUp[index][1],:] = [0,255,0]
            if coordUp[index][0]-1 >= 0 and coordUp[index][0]-1 < USimageP.shape[0]:
                USimageP[coordUp[index][0]-1, coordUp[index][1],:] = [0,255,0]
        if coordDeep[index][0] >= 0 and coordDeep[index][0] < USimageP.shape[0]:
            USimageP[coordDeep[index][0], coordDeep[index][1], :] = [0, 255, 0]
            if coordDeep[index][0]+1 >= 0 and coordDeep[index][0]+1 < USimageP.shape[0]:
                USimageP[coordDeep[index][0]+1, coordDeep[index][1], :] = [0,255,0]
            if coordDeep[index][0]-1 >= 0 and coordDeep[index][0]-1 < USimageP.shape[0]:
                USimageP[coordDeep[index][0]-1, coordDeep[index][1], :] = [0,255,0]  

    #snippets
    for f in range(len(fasc)):
        for g in range(fasc[f].shape[0]):
            USimageP[fasc[f][g,0], fasc[f][g,1], :] = [255,255,255]
    
    #fascicles
    PAu = []
    PAd = []
    outPAu = []
    outPAd = []
    fasc_in = []
    fasc_out = []
    for n1 in range(len(splines_fasc)):
        newy = np.arange(0, USimageP.shape[1], 1)
        newx = np.int32(splines_fasc[n1](newy))
        coord = np.vstack((newx, newy)).T
        if muscle_fascicles[n1][8] == 'in image':
            for n2 in range(coord.shape[0]):
                if coord[n2][0]>=0 and coord[n2][0]<USimageP.shape[0]:
                    USimageP[coord[n2][0], coord[n2][1], :] = [0,255,0]
            PAu.append(muscle_fascicles[n1][4])
            PAd.append(muscle_fascicles[n1][7])
            fasc_in.append(muscle_fascicles[n1][9])
            
        else:
            outPAu.append(muscle_fascicles[n1][4])
            outPAd.append(muscle_fascicles[n1][7])
            fasc_out.append(muscle_fascicles[n1][9])
            for n2 in range(coord.shape[0]):
                if coord[n2][0]>=0 and coord[n2][0]<USimageP.shape[0]:
                    USimageP[coord[n2][0], coord[n2][1], :] = [0,0,255]
#                if coord[n2][0]+1>=0 and coord[n2][0]+1<USimageP.shape[0]:
#                    USimageP[coord[n2][0]+1, coord[n2][1], :] = [0,0,255]
#                if coord[n2][0]-1>=0 and coord[n2][0]-1<USimageP.shape[0]:  
#                    USimageP[coord[n2][0]-1, coord[n2][1], :] = [0,0,255]      
        
    cv2.imshow('full image', USimageP)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()

    #FEATURES
    import matplotlib.pyplot as plt
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharex = False, sharey = False, figsize =(15,15))
    
    #muscle thickness
    ax1.plot(abscissa, thickness, 'b+', markersize = 3)
    ax1.set_title('Muscle thickness')
    ax1.set_ylabel('thickness (mm)')
    ax1.set_xlabel('muscle axis (mm)')
    ax1.set_ybound(0,30)
    
    #pennation angle
    ax2.plot(PAu, PAd, 'g+', markersize = 5, label = 'inside image')
    ax2.plot(outPAu, outPAd, 'r+', markersize = 5, label = 'extrapolated')
    plt.legend(loc = 'upper left', prop={'size': 6})
    ax2.set_title('Pennation angles')
    ax2.set_xlabel('Angle (degrees) / upper apo')
    ax2.set_ylabel('Angle (degrees) / deep apo')
    ax2.set_xbound(0,30)
    ax2.set_ybound(0,30)

    #fascicles length
    ax3.plot(fasc_in, [1] * len(fasc_in), 'g+', markersize = 5, label = 'inside image')
    ax3.plot(fasc_out, [1] * len(fasc_out), 'r+', markersize = 5, label = 'extrapolated')
    plt.legend(loc = 'upper left', prop={'size': 6})
    ax3.set_title('Fascicles length')
    ax3.set_xlabel('length (mm)')
    
    plt.show()
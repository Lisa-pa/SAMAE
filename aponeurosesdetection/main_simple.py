from calibration.calib import autoCalibration
from preprocessing.cropping import autocropping
from preprocessing.preprocess import preprocessingApo
import aponeuroseslocation as apoL
import aponeuroses_contours as apoC
import MUFeaM
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
RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/31_romain/architecture/Romain_jamon_20181008_084409_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/fam_1/architecture/Julien_jamon_20180720_170728_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/34_nicolas/architecture/Nicolas_post_20181210_105644_image_bfs.jpg', -1)
#RGBimage = data.simpleimg()
#################################################
cv2.imshow('Image to process', RGBimage)
process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
    Close all windows after clicking yes or no.', default = 'yes', icon='question')
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()

if process == True:
    #Calibrate the image
    # calibX, calibY are the calibration factors in X and Y directions
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
            entry5 = input("Please enter value in millimeters representing the optional additional cropping from top of the image:")
            entry6 = input("Please enter value in millemeters representing the optional additional cropping from bottom of the image:")
            THRESH1 = int(entry1)
            THRESH2 = int(entry2)
            THRESH3 = int(entry3)
            THRESH4 = int(entry4)
            CROP1 = int(entry5)
            CROP2 = int(entry6)
            print(f'You entered the following thresholds: {THRESH1, THRESH2, THRESH3, THRESH4}')
            if THRESH1>255 or THRESH1<0 or THRESH2>255 or THRESH2<0 or THRESH3>255 or THRESH3<0 or\
                THRESH4>255 or THRESH4<0:
                raise ValueError('All thresholds must be integers between 0 and 255')
            USimage = autocropping(RGBimage, THRESH1, THRESH2, THRESH3, THRESH4, calibY, CROP1, CROP2)
#    cv2.imwrite('C:/Users/Lisa Paillard/Desktop/julien728cropped.jpg', USimage)
    #################################################

    #Locate aponeuroses and find linear approximation of aponeuroses

    # USimage_pp:   cropped ulstrasound image that underwent pre-processing
    # paramUp:      parameters (a,b) for the line equation x = ay + b 
    #               that linearly approximates superficial aponeurosis. 
    #               x = rows, y = columns
    # paramLow:     parameters for the line equation of deep aponeurosis
    # apoUp:        (l1,l2) are the two lines in between which the 
    #               upper aponeurosis spotted.
    # apoUp:        lines in between which deep aponeurosis was spotted
    USimage_pp = preprocessingApo(USimage, 'localmean', 0, 41) #pre_processing
    paramUp, paramLow, apoUp, apoLow = apoL.twoApoLocation(USimage_pp, thresh = None)
    
    '''USimage[apoUp[0], :,:] = [0,0,255]
    USimage[apoUp[1], :,:] = [0,0,255]
    USimage[apoLow[0], :,:] = [0,255,0]
    USimage[apoLow[1], :,:] = [0,255,0]
    cv2.imshow('interpolated aponeuroses', USimage)
    cv2.imshow('pp image', USimage_pp)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()'''
    
    #sub-images containing superficial and deep aponeurosis each
    upperApo = np.copy(USimage[apoUp[0]:apoUp[1],:])
    lowerApo = np.copy(USimage[apoLow[0]:apoLow[1],:])
    
    #Pre-processed sub-images of superficial and deep aponeuroses
    upperApo_pp = np.copy(USimage_pp[apoUp[0]:apoUp[1],:])
    lowerApo_pp = np.copy(USimage_pp[apoLow[0]:apoLow[1],:])
    
    #################################################
    
    # Get exact contour of each aponeurosis
   
    ini_upApo = apoC.initiateContour(upperApo_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1] - apoUp[0], 8])
    contourUp, nUp = apoC.activeContour(upperApo_pp,ini_upApo,0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
    print('Upper aponeurosis contour found in ', nUp, ' steps')

    if np.min(contourUp)>0: #it means active contour model failed
        #try a second time with another initial contour
        ini_upApo = apoC.initiateContour(upperApo_pp, typeC = 'quadrangle_param', param = [paramUp[0], paramUp[1] - apoUp[0], 40])
        contourUp, nUp = apoC.activeContour(upperApo_pp, ini_upApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
        print('Upper aponeurosis contour found in ', nUp, ' steps')

        #if min(contourUp) still positive
        #use linear approximation because it means active contour model failed again
        if np.min(contourUp)>0:
            yUp = np.arange(0,USimage.shape[1],1)
            xUp = np.int32(yUp*paramUp[0]+paramUp[1])
            coordUp = np.vstack((xUp, yUp)).T

    if np.min(contourUp)<=0: #ask for validation of the contour
        contourUpimage, contourPointsUp = apoC.extractContour(contourUp, upperApo, offSetX = apoUp[0], offSetY = 0)
        cv2.imshow('Upper aponeurosis contour', contourUpimage)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        
        if valid == True: #B-spline to approximate aponeurosis if contour suits
            approxUp, coordUp = apoC.approximate(contourPointsUp, 'upper', upperApo_pp, d = 1)
        elif valid == False: #use linear approximation from radon transform
            yUp = np.arange(0,USimage.shape[1],1)
            xUp = np.int32(yUp*paramUp[0]+paramUp[1])
            coordUp = np.vstack((xUp, yUp)).T


    ini_lowApo = apoC.initiateContour(lowerApo_pp, typeC = 'quadrangle_param', param = [paramLow[0], paramLow[1] - apoLow[0], 8])
    contourLow, nLow = apoC.activeContour(lowerApo_pp, ini_lowApo, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
    print('Lower aponeurosis contour found in ', nLow, ' steps')
    
    if np.min(contourLow)>0: #it means active contour model failed
        #try a second time with another initial contour
        ini_lowApo = apoC.initiateContour(lowerApo_pp, typeC = 'quadrangle_param', param = [paramLow[0], paramLow[1] - apoLow[0], 40])
        contourLow, nLow = apoC.activeContour(lowerApo_pp, ini_lowApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
        print('Lower aponeurosis contour found in ', nUp, ' steps')

        #if min(contourLow) still positive
        #use linear approximation because it means active contour model failed again
        if np.min(contourLow)>0:
            yLow = np.arange(0,USimage.shape[1],1)
            xLow = np.int32(yLow*paramLow[0]+paramLow[1])
            coordLow = np.vstack((xLow, yLow)).T    

    if np.min(contourLow)<=0: #ask for validation of the contour
        contourLowimage, contourPointsLow = apoC.extractContour(contourLow, lowerApo, offSetX = apoLow[0], offSetY = 0)
        cv2.imshow('Lower aponeurosis contour', contourLowimage)
        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()    
        
        if valid == True: #B-spline to approximate aponeurosis if contour suits
            approxLow, coordLow = apoC.approximate(contourPointsLow, 'lower', lowerApo_pp, d = 1)
        elif valid == False:#else linear approximation from Radon transform is used
            yLow = np.arange(0,USimage.shape[1],1)
            xLow = np.int32(yLow*paramLow[0]+paramLow[1])
            coordLow = np.vstack((xLow, yLow)).T
    
    
    #Muscle thickness calculation
    absc, thickness, spline_thickness = MUFeaM.muscleThickness(USimage, coordUp, coordLow, 0, USimage.shape[1]-1, calibX, calibY)
    
    #Crop to focus in between aponeuroses
    #ROI : region of interest
    crop1 = np.amax(coordUp[:,0]) + 1 
    crop2 = np.amin(coordLow[:,0])
    ROI = np.copy(USimage[crop1:crop2, :, :])
    
    #Enhance tube-like structures with MVEF method - Frangi - 
    #let's consider that fascicle diameter is between 0.3 mm and 0.5 mm,
    #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
    sca = np.arange(round(0.3/calibX), round(0.5/calibX), 0.5)
    print('MVEF running')
    MVEF_image = FaDe.MVEF_2D(255-ROI, sca, [0.5, 0])
    #
    #threshold
    threshMVEF_percent = 85
    threshMVEF = np.percentile(MVEF_image, threshMVEF_percent)
    MVEF_image2 = cv2.threshold(MVEF_image, threshMVEF, 255, cv2.THRESH_BINARY)[1]

    #locate snippets and filter them to locate muscle fascicles
    snippets, snip_lines = FaDe.locateSnippets(MVEF_image2, calibX, calibY,\
                                               minLength = 4, rangeAngles = [6,45],\
                                               percentageAlign = 0.70, offSetX = crop1)
    fasc = FaDe.combineSnippets(USimage, snippets, snip_lines, thresh_alignment = 10, thresh_length = 200)

    #transform the snippets, which are in fact the contours of the fascicles,
    #into lines by taking the mean of the contour, that is the mean on each
    #column of the contour. When branches exist, the mean is not computed (the
    #region is ignored)
    averages = [] 
    for i in range(len(fasc)):
        averages.append(FaDe.contourAverage(fasc[i]))
    
    #interpolations to get fascicles' curve
    interpolated_fasc = FaDe.approximateFascicle(USimage, averages, d = 1)
    splines = [interpolated_fasc[i][0] for i in range(len(interpolated_fasc))]
    
    #intersections of fascicles with aponeuroses
    intersecL = MUFeaM.findIntersections(approxLow, splines, USimage, typeI = 'simple')
    intersecU = MUFeaM.findIntersections(approxUp, splines, USimage, typeI = 'simple')
    
    #Pennation angles calculation
    PA_up = MUFeaM.pennationAngles(approxUp, splines, intersecU, calibX, calibY, USimage)
    PA_low = MUFeaM.pennationAngles(approxLow, splines, intersecL, calibX, calibY, USimage)
    print('angles up', PA_up)
    print('angles low', PA_low)
    
    #Fascicles length calculation
    fasc_length = MUFeaM.fasciclesLength(splines, intersecU, intersecL, calibX, calibY)
    print(fasc_length)
    
    
    ################
    #Visualization:

    Zerosvertic = np.uint8(np.zeros(USimage.shape))
    ImF = cv2.hconcat([Zerosvertic,Zerosvertic, USimage, Zerosvertic, Zerosvertic])

    #aponeuroses
    newy = np.arange(-2*USimage.shape[1], 3*USimage.shape[1], 1)
    newx_u = np.int32(approxUp(newy))
    newx_l = np.int32(approxLow(newy))
    newy = newy + 2*USimage.shape[1]
    coordUp = np.vstack((newx_u, newy)).T
    coordLow = np.vstack((newx_l, newy)).T
    
    for index in range(coordUp.shape[0]):
        if coordUp[index][0] >= 0 and coordUp[index][0] < ImF.shape[0]:
            ImF[coordUp[index][0], coordUp[index][1],:] = [255,0,0]
            if coordUp[index][0] +1 >= 0 and coordUp[index][0]+1 < ImF.shape[0]:
                ImF[coordUp[index][0]+1, coordUp[index][1],:] = [255,0,0]
            if coordUp[index][0]-1 >= 0 and coordUp[index][0]-1 < ImF.shape[0]:
                ImF[coordUp[index][0]-1, coordUp[index][1],:] = [255,0,0]
        if coordLow[index][0] >= 0 and coordLow[index][0] < ImF.shape[0]:
            ImF[coordLow[index][0], coordLow[index][1], :] = [255,0,0]
            if coordLow[index][0]+1 >= 0 and coordLow[index][0]+1 < ImF.shape[0]:
                ImF[coordLow[index][0]+1, coordLow[index][1], :] = [255,0,0]
            if coordLow[index][0]-1 >= 0 and coordLow[index][0]-1 < ImF.shape[0]:
                ImF[coordLow[index][0]-1, coordLow[index][1], :] = [255,0,0]  

    #fascicles
    for f in range(len(fasc)):
        for g in range(fasc[f].shape[0]):
            ImF[fasc[f][g,0], fasc[f][g,1] + 2*USimage.shape[1], :] = [255,255,255]
    
    for a in range(len(interpolated_fasc)):
#        coord = interpolated_fasc[a][1]
        newy = np.arange(-2*USimage.shape[1], 3*USimage.shape[1], 1)
        newx = np.int32(interpolated_fasc[a][0](newy))
        newy = newy + 2*USimage.shape[1]
        coord = np.vstack((newx, newy)).T
        for b in range(coord.shape[0]):
            if coord[b][0]>=0 and coord[b][0]<ImF.shape[0]:
                ImF[coord[b][0], coord[b][1], :] = [0,255,0]
#                if coord[b][0]+1>=0 and coord[b][0]+1<ImF.shape[0]:
#                    ImF[coord[b][0]+1, coord[b][1], :] = [0,255,0]
#                if coord[b][0]-1>=0 and coord[b][0]-1<ImF.shape[0]:  
#                    ImF[coord[b][0]-1, coord[b][1], :] = [0,255,0]
    
       
    #intersection points
    for i1 in range(len(splines)):
        xi = int(intersecL[i1][0])
        yi = int(intersecL[i1][1]) + 2*USimage.shape[1]
    
        if xi>=3 and xi<ImF.shape[0]-3\
        and yi>=3 and yi<ImF.shape[1]-3:
            ImF[xi,yi,:] = [255,255,50]
            ImF[xi+1,yi,:] = [255,255,50]
            ImF[xi-1,yi,:] = [255,255,50]
            ImF[xi,yi+1,:] = [255,255,50]
            ImF[xi,yi-1,:] = [255,255,50]
            ImF[xi+2,yi,:] = [255,255,50]
            ImF[xi-2,yi,:] = [255,255,50]
            ImF[xi,yi+2,:] = [255,255,50]
            ImF[xi,yi-2,:] = [255,255,50]
            ImF[xi+3,yi,:] = [255,255,50]
            ImF[xi-3,yi,:] = [255,255,50]
            ImF[xi,yi+3,:] = [255,255,50]
            ImF[xi,yi-3,:] = [255,255,50]
            ImF[xi+1,yi-3,:] = [255,255,50]
            ImF[xi+1,yi-2,:] = [255,255,50]
            ImF[xi+1,yi-1,:] = [255,255,50]
            ImF[xi+1,yi+3,:] = [255,255,50]
            ImF[xi+1,yi+2,:] = [255,255,50]
            ImF[xi+1,yi+1,:] = [255,255,50] 
            ImF[xi-1,yi-3,:] = [255,255,50]
            ImF[xi-1,yi-2,:] = [255,255,50]
            ImF[xi-1,yi-1,:] = [255,255,50]
            ImF[xi-1,yi+3,:] = [255,255,50]
            ImF[xi-1,yi+2,:] = [255,255,50]
            ImF[xi-1,yi+1,:] = [255,255,50]
            ImF[xi-3,yi+1,:] = [255,255,50] 
            ImF[xi-2,yi+1,:] = [255,255,50] 
            ImF[xi+3,yi+1,:] = [255,255,50] 
            ImF[xi+2,yi+1,:] = [255,255,50] 
            ImF[xi-3,yi-1,:] = [255,255,50]             
            ImF[xi-2,yi-1,:] = [255,255,50]             
            ImF[xi+3,yi-1,:] = [255,255,50]             
            ImF[xi+2,yi-1,:] = [255,255,50]
            
        xi = int(intersecU[i1][0])
        yi = int(intersecU[i1][1]) + 2*USimage.shape[1]
    
        if xi>=3 and xi<ImF.shape[0]-3\
        and yi>=3 and yi<ImF.shape[1]-3:
            ImF[xi,yi,:] = [255,255,50]
            ImF[xi+1,yi,:] = [255,255,50]
            ImF[xi-1,yi,:] = [255,255,50]
            ImF[xi,yi+1,:] = [255,255,50]
            ImF[xi,yi-1,:] = [255,255,50]
            ImF[xi+2,yi,:] = [255,255,50]
            ImF[xi-2,yi,:] = [255,255,50]
            ImF[xi,yi+2,:] = [255,255,50]
            ImF[xi,yi-2,:] = [255,255,50]
            ImF[xi+3,yi,:] = [255,255,50]
            ImF[xi-3,yi,:] = [255,255,50]
            ImF[xi,yi+3,:] = [255,255,50]
            ImF[xi,yi-3,:] = [255,255,50]
            ImF[xi+1,yi-3,:] = [255,255,50]
            ImF[xi+1,yi-2,:] = [255,255,50]
            ImF[xi+1,yi-1,:] = [255,255,50]
            ImF[xi+1,yi+3,:] = [255,255,50]
            ImF[xi+1,yi+2,:] = [255,255,50]
            ImF[xi+1,yi+1,:] = [255,255,50] 
            ImF[xi-1,yi-3,:] = [255,255,50]
            ImF[xi-1,yi-2,:] = [255,255,50]
            ImF[xi-1,yi-1,:] = [255,255,50]
            ImF[xi-1,yi+3,:] = [255,255,50]
            ImF[xi-1,yi+2,:] = [255,255,50]
            ImF[xi-1,yi+1,:] = [255,255,50]
            ImF[xi-3,yi+1,:] = [255,255,50] 
            ImF[xi-2,yi+1,:] = [255,255,50] 
            ImF[xi+3,yi+1,:] = [255,255,50] 
            ImF[xi+2,yi+1,:] = [255,255,50] 
            ImF[xi-3,yi-1,:] = [255,255,50]             
            ImF[xi-2,yi-1,:] = [255,255,50]             
            ImF[xi+3,yi-1,:] = [255,255,50]             
            ImF[xi+2,yi-1,:] = [255,255,50]             
    
    ImF = cv2.resize(src = ImF, dsize = (int(ImF.shape[1]*0.4), int(ImF.shape[0]*0.4)), interpolation = cv2.INTER_CUBIC)
    cv2.imshow('Final image', ImF)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()

    #FEATURES
    import matplotlib.pyplot as plt
    listAbs = [1 for ab in range(len(splines))]
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharex = False, sharey = False, figsize =(15,15))
    
    #muscle thickness
    ax1.plot(absc, thickness, 'b+', markersize = 3)
    ax1.set_title('Muscle thickness')
    ax1.set_ylabel('thickness (mm)')
    ax1.set_xlabel('muscle axis (mm)')
    ax1.set_ybound(0,30)
    
    #pennation angle
    ax2.plot(PA_up, PA_low, 'r+', markersize = 5)
    ax2.set_title('Pennation angles')
    ax2.set_xlabel('Angle (degrees) / upper apo')
    ax2.set_ylabel('Angle (degrees) / deep apo')
    ax2.set_xbound(0,30)
    ax2.set_ybound(0,30)

    #fascicles length
    ax3.plot(fasc_length, listAbs, 'k+', markersize = 5)
    ax3.set_title('Fascicles length')
    ax3.set_xlabel('length (mm)')
    
    plt.show()
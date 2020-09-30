from calibration.calib import autoCalibration
from preprocessing.cropping import autocropping
from preprocessing.preprocess import preprocessingApo
import apoLoc as apoL
import apoCont as apoC
import MUFeaM
import FaDe as FaDe

import cv2
import numpy as np
import tkinter.messagebox as tkbox


########################PROCESSING OF SIMPLE US IMAGES#########################

#Open the image
#import aponeurosesdetection.data as apodat

#! RGBimage = apodat.simpleimg()
#RGBimage = cv2.imread('C:/Users/Antonio/Desktop/AponeurosesDetection/aponeurosesdetection/data/simple_echo.jpg', -1)

#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/01_Kevin/POST/architecture/Kevin_jamon_20180927_142814_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/31_romain/architecture/Romain_jamon_20181008_084433_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/fam_1/architecture/Julien_jamon_20180720_170728_image_bfs.jpg', -1)
#RGBimage = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/data/34_nicolas/architecture/Nicolas_post_20181210_105644_image_bfs.jpg', -1)
#RGBimage = data.simpleimg()



def simpleprocessing(path_to_img):

    RGBimage = cv2.imread(path_to_img, -1)
    
    #Validate the image processing start
    cv2.imshow('Image to process', RGBimage)
    process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
        Close all windows after clicking yes or no.', default = 'yes', icon='question')
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    
    if process == True:
        

        #Calibrate the image
        # calibX, calibY are the calibration factors in X and Y directions
        calibX, calibY = autoCalibration(RGBimage)
        if calibX > 2 * calibY or calibY > 2 * calibX:
            calibX = min(calibX, calibY)
            calibY = min(calibX, calibY)
                    
        #Crop the image to keep essential US data
        USimage, l1, l2, c1, c2 = autocropping(RGBimage, 10., 15., 12., 25., calibY, additionalCrop1 = 2., additionalCrop2 = 6.)

        #--------ask for validation; if cropping does not suit the user, ask the user for thresholds
        ok = False
        counter = 1
        while ok == False and counter<=5:
            if USimage.size>0:
                cv2.imshow('Cropped US image', USimage)
                ok = tkbox.askyesno('Need user validation', 'Do you validate the cropping ? If no, we will ask you new thresholds in the command prompt. After clicking yes or no, please close the image window to continue.', default = 'yes', icon='question')
                cv2.waitKey(0) & 0xFF
                cv2.destroyAllWindows()
    
            if ok == False:
                counter = counter + 1
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
                USimage, l1, l2, c1, c2 = autocropping(RGBimage, THRESH1, THRESH2, THRESH3, THRESH4, calibY, CROP1, CROP2)

        if counter > 5:
            archi_auto = dict()
            archi_auto['crop'] = {'lines': 'error', 'columns': 'error'}
            return archi_auto
        
        cv2.imwrite(path_to_img[:-8]+'_cropped.jpg', USimage)    

        #Locate aponeuroses and find linear approximation of aponeuroses
    
        # USimage_pp:   cropped ulstrasound image that underwent pre-processing
        # paramSup:      parameters (a,b) for the line equation x = ay + b 
        #               that linearly approximates superficial aponeurosis. 
        #               x = rows, y = columns
        # paramInf:     parameters for the line equation of deep aponeurosis
        # locSup:        (l1,l2) are the two lines in between which the 
        #               upper aponeurosis spotted.
        # locInf:        lines in between which deep aponeurosis was spotted
        USimage_pp = preprocessingApo(I = USimage, typeI = 'simple', mode = 'localmean', margin = 0, sizeContrast = 41) #pre_processing
        paramSup, paramInf, locSup, locInf = apoL.twoApoLocation(USimage_pp, angle1 = 80, angle2 = 100, thresh = None, calibV = calibX)
        cv2.imwrite(path_to_img[:-8]+'_preprocessed.jpg', USimage_pp)        
        
        
        if paramSup[0] == 'error':
            archi_auto = dict()
            archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
            archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
            archi_auto['aposup'] = {'coords':'error'}
            archi_auto['apoinf'] = {'coords':'error'}
            return archi_auto
        
        '''
        USimage[locSup[0], :,:] = [0,0,255]
        USimage[locSup[1], :,:] = [0,0,255]
        USimage[locInf[0], :,:] = [0,255,0]
        USimage[locInf[1], :,:] = [0,255,0]
        cv2.imshow('interpolated aponeuroses', USimage)
        cv2.imshow('pp image', USimage_pp)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        '''
        
        #sub-images containing superficial and deep aponeurosis each
        SupApo = np.copy(USimage[locSup[0]:locSup[1],:])
        InfApo = np.copy(USimage[locInf[0]:locInf[1],:])
        
        #Pre-processed sub-images of superficial and deep aponeuroses
        SupApo_pp = np.copy(USimage_pp[locSup[0]:locSup[1],:])
        InfApo_pp = np.copy(USimage_pp[locInf[0]:locInf[1],:])

        # Get exact contour of each aponeurosis
        iniSupApo = apoC.initiateContour(SupApo_pp, typeC = 'quadrangle_param', param = [paramSup[0], paramSup[1] - locSup[0], 8])
        contourSup, nSup = apoC.activeContour(SupApo_pp,iniSupApo,0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
        print('Upper aponeurosis contour found in ', nSup, ' steps')
    
        if np.min(contourSup)>0: #it means active contour model failed
            #try a second time with another initial contour
            iniSupApo = apoC.initiateContour(SupApo_pp, typeC = 'quadrangle_param', param = [paramSup[0], paramSup[1] - locSup[0], 40])
            contourSup, nSup = apoC.activeContour(SupApo_pp, iniSupApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
            print('Upper aponeurosis contour found in ', nSup, ' steps')
            #if min(contourSup) still positive
            #use linear approximation because it means active contour model failed again
            type_approx_UA = 'linear'
    
        if np.min(contourSup)<=0: #ask for validation of the contour
            contourSup_image, contourSup_points = apoC.extractContour(contourSup, SupApo, offSetX = locSup[0], offSetY = 0)
            cv2.imshow('Upper aponeurosis contour', contourSup_image)
            valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
            cv2.waitKey(0) & 0xFF
            cv2.destroyAllWindows()
            
            if valid == False:
                iniSupApo = apoC.initiateContour(SupApo_pp, typeC = 'quadrangle_param', param = [0, int(SupApo.shape[0]/2), 40])
                contourSup, nSup = apoC.activeContour(SupApo_pp, iniSupApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
                print('Upper aponeurosis contour found in ', nSup, ' steps')
                if np.min(contourSup)<=0:
                    contourSup_image, contourSup_points = apoC.extractContour(contourSup, SupApo, offSetX = locSup[0], offSetY = 0)
                    cv2.imshow('Upper aponeurosis contour', contourSup_image)
                    valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
                    cv2.waitKey(0) & 0xFF
                    cv2.destroyAllWindows()                
                    
            if valid == True:   #B-spline to approximate aponeurosis if contour suits
                                #replace paramSup coefficients by the spline
                type_approx_UA = 'spline'
                paramSup = apoC.approximateApo(p = contourSup_points, apoType = 'upper', I = SupApo_pp, typeapprox = 'polyfit', d = 1)
            elif valid == False: #use linear approximation from radon transform
                type_approx_UA = 'linear'
    
    
    
        ini_InfApo = apoC.initiateContour(InfApo_pp, typeC = 'quadrangle_param', param = [paramInf[0], paramInf[1] - locInf[0], 8])
        contourInf, nInf = apoC.activeContour(InfApo_pp, ini_InfApo, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
        print('Lower aponeurosis contour found in ', nInf, ' steps')
        
        if np.min(contourInf)>0: #it means active contour model failed
            #try a second time with another initial contour
            ini_InfApo = apoC.initiateContour(InfApo_pp, typeC = 'quadrangle_param', param = [paramInf[0], paramInf[1] - locInf[0], 40])
            contourInf, nInf = apoC.activeContour(InfApo_pp, ini_InfApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
            print('Lower aponeurosis contour found in ', nInf, ' steps')
            #if min(contourInf) still positive
            #use linear approximation because it means active contour model failed again
            type_approx_LA = 'linear'
    
        if np.min(contourInf)<=0: #ask for validation of the contour
            contourInf_image, contourInf_Points = apoC.extractContour(contourInf, InfApo, offSetX = locInf[0], offSetY = 0)
            cv2.imshow('Lower aponeurosis contour', contourInf_image)
            valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
            cv2.waitKey(0) & 0xFF
            cv2.destroyAllWindows()    
            if valid == False:
                ini_InfApo = apoC.initiateContour(InfApo_pp, typeC = 'quadrangle_param', param = [0, int(InfApo.shape[0]/2), 40])
                contourInf, nInf = apoC.activeContour(InfApo_pp, ini_InfApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
                print('Lower aponeurosis contour found in ', nInf, ' steps')
                
                if np.min(contourInf)<=0:
                    contourInf_image, contourInf_Points = apoC.extractContour(contourInf, InfApo, offSetX = locInf[0], offSetY = 0)
                    cv2.imshow('Lower aponeurosis contour', contourInf_image)
                    valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
                    cv2.waitKey(0) & 0xFF
                    cv2.destroyAllWindows()    

            if valid == True:   #B-spline to approximate aponeurosis if contour suits
                                #replace paramSup coefficients by the spline
                type_approx_LA = 'spline'
                paramInf = apoC.approximateApo(p = contourInf_Points, apoType = 'lower', I = InfApo_pp, typeapprox = 'polyfit', d = 1)
            elif valid == False:#else linear approximation from Radon transform is used
                type_approx_LA = 'linear'
        
    
        #calculate coordinates of aponeuroses
        #if aponeusosis linear, param is transformed from a list of parameters into a spline
        coordSup, paramSup = MUFeaM.pointsCoordinates(typeA = type_approx_UA, param = paramSup, interval = [0, USimage.shape[1]])
        coordInf, paramInf = MUFeaM.pointsCoordinates(typeA = type_approx_LA, param = paramInf, interval = [0, USimage.shape[1]])
    
       
        #Muscle thickness calculation
        absc, thickness, spline_thickness = MUFeaM.muscleThickness(spl1 = paramSup, spl2 = paramInf, start = 0, end = USimage.shape[1]-1, calibV = calibX, calibH = calibY)
            
            
        #ROI : region of interest in between aponeuroses
        crop1 = np.amax(coordSup[:,0]) + 30
        crop2 = np.amin(coordInf[:,0]) - 30
        ROI = np.copy(USimage[crop1:crop2, :, :])
        cv2.imwrite(path_to_img[:-8]+'_ROI.jpg', ROI)
                
        
        #Enhance tube-like structures with MVEF method - Frangi - 
        #let's consider that fascicle diameter is between 0.3 mm and 0.5 mm,
        #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
        sca = np.arange(round(0.3/calibX), round(0.5/calibX), 0.5)
        print('MVEF running')
        MVEF_image = FaDe.MVEF_2D(255-ROI, sca, [0.5, 0])
        cv2.imwrite(path_to_img[:-8]+'_MVEF.jpg', MVEF_image)

        #
        #threshold
        threshMVEF_percent = 85
        threshMVEF = np.percentile(MVEF_image, threshMVEF_percent)
        MVEF_image2 = cv2.threshold(MVEF_image, threshMVEF, 255, cv2.THRESH_BINARY)[1]
        cv2.imwrite(path_to_img[:-8]+'_MVEF2.jpg', MVEF_image2)
    
        #locate snippets and filter them to locate muscle fascicles
        snippets, snip_lines = FaDe.locateSnippets(MVEF_image2, calibX, calibY,\
                                                   minLength = 4,\
                                                   offSetX = crop1, im = USimage)
        
        if snippets == 'error':
            archi_auto = dict()
            archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
            archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
            archi_auto['aposup'] = {'coords':coordSup}
            archi_auto['apoinf'] = {'coords':coordInf}
            archi_auto['MT'] = {'coords': np.vstack((absc, thickness)).T, 'columns interval': [0, USimage.shape[1]-1]}
            return archi_auto 
        
        
        fasc, fasc2 = FaDe.combineSnippets(USimage, snippets, snip_lines, min_nb_sn = 2, thresh_alignment = 20)
        print('fascicles', len(fasc))
        print ('fascicles2', len(fasc2))
    
    
        #transform the snippets, which are in fact the contours of the fascicles,
        #into lines by taking the mean line inside the contour
        #When branches exist, the mean is not computed (the region is ignored)
        averages = [] 
        for i in range(len(fasc)):
            averages.append(FaDe.contourAverage(fasc[i]))

        averages2 = [] 
        for i in range(len(fasc2)):
            averages2.append(FaDe.contourAverage(fasc2[i]))        
        
        
        #interpolations to get fascicles' curve
        splines_fasc = FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages, d = 1)
        splines_fasc = splines_fasc + FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages2, d = 1)
        
        #intersections of fascicles with aponeuroses (in pixels)
        intersecL, intersecU, splines_fasc = MUFeaM.findIntersections(spl_inf = paramInf, spl_sup = paramSup, listSpl = splines_fasc, start = USimage.shape[1]/2, insertion = 200/calibY)
        
        #Pennation angles calculation (in degrees)
        PASup = MUFeaM.pennationAngles(paramSup, splines_fasc, intersecU, calibX, calibY)
        PA_Inf = MUFeaM.pennationAngles(paramInf, splines_fasc, intersecL, calibX, calibY)
        
        #Fascicles length calculation (in millimeters)
        fasc_length = MUFeaM.fasciclesLength(splines_fasc, intersecU, intersecL, calibX, calibY)
        
        #fascicle location in the original RGB image
        loc_fasc = MUFeaM.locateFasc(intersecL, [-l1,-c1], calibY)
        
        #Dict containing info per image
        archi_auto = dict()
        archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
        archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
        archi_auto['aposup'] = {'coords':coordSup}
        archi_auto['apoinf'] = {'coords':coordInf}
        archi_auto['MT'] = {'coords': np.vstack((absc, thickness)).T, 'columns interval': [0, USimage.shape[1]-1]}
        
        for index in range(len(splines_fasc)):
            archi_auto['fsc_' + str(index+1)] = {}
            typeIU = 'out of image'
            typeIL = 'out of image'
            typeFL = 'out of image'
            if intersecU[index][0] >= 0 and\
                intersecU[index][0] < USimage.shape[0] and\
                intersecU[index][1] >= 0 and\
                intersecU[index][1] < USimage.shape[1]:
                typeIU = 'in image'
            if intersecL[index][0] >= 0 and\
                intersecL[index][0] < USimage.shape[0] and\
                intersecL[index][1] >= 0 and\
                intersecL[index][1] < USimage.shape[1]:
                typeIL = 'in image'
            if typeIU == 'in image' and typeIL == 'in image':
                typeFL = 'in image'           
            archi_auto['fsc_' + str(index+1)] = {'dist from (0,0) of RGB image, in mm': loc_fasc[index]}
            archi_auto['fsc_' + str(index+1)]['PAsup'] = {'in/out of the image': typeIU,
                                                          'value in degree' : PASup[index],
                                                          'intersection with apo': intersecU[index]}

            archi_auto['fsc_' + str(index+1)]['PAinf'] = {'in/out of the image': typeIL,
                                                          'value in degree' : PA_Inf[index],
                                                          'intersection with apo': intersecL[index]}
            
            archi_auto['fsc_' + str(index+1)]['FL'] = {'in/out of the image':typeFL,
                                                       'length in mm': fasc_length[index]}   
    
    
    
    
    
    
    
    
    
    
    
        ################
        #Visualization:
        #snippets
        couleurs = [[255,0,0], [0,255,0], [0,0,255], [255,255,0],[255,0,255], [0,255,255],\
                    [100,200,0],[100,200,100], [50,200,0],[50,100,50], [255,100,0],\
                    [120,120,255], [255,80,80],[0,100,200], [0,100,80], [255,255,255],\
                    [120,120,120], [50,100,150],[100,50,150], [150,100,50], [50,150,100],
                    [100,150,50],[150,50,100],[12,75,255],[40,140,40]]
        for f in range(len(fasc)):
            for g in range(fasc[f].shape[0]):
                USimage[fasc[f][g,0], fasc[f][g,1], :] = couleurs[f]
        for f in range(len(fasc2)):
            for g in range(fasc2[f].shape[0]):
                USimage[fasc2[f][g,0], fasc2[f][g,1], :] = couleurs[-f]
                
        Zerosvertic = np.uint8(np.zeros(USimage.shape))
        ImF = cv2.hconcat([Zerosvertic,Zerosvertic, USimage, Zerosvertic, Zerosvertic])
    
        #aponeuroses
        newy = np.arange(-2*USimage.shape[1], 3*USimage.shape[1], 1)
        newx_u = np.int32(paramSup(newy))
        newx_l = np.int32(paramInf(newy))
        newy = newy + 2*USimage.shape[1]
        coordSup = np.vstack((newx_u, newy)).T
        coordInf = np.vstack((newx_l, newy)).T
        
        for index in range(coordSup.shape[0]):
            if coordSup[index][0] >= 0 and coordSup[index][0] < ImF.shape[0]:
                ImF[coordSup[index][0], coordSup[index][1],:] = [255,0,0]
                if coordSup[index][0] +1 >= 0 and coordSup[index][0]+1 < ImF.shape[0]:
                    ImF[coordSup[index][0]+1, coordSup[index][1],:] = [255,0,0]
                if coordSup[index][0]-1 >= 0 and coordSup[index][0]-1 < ImF.shape[0]:
                    ImF[coordSup[index][0]-1, coordSup[index][1],:] = [255,0,0]
            if coordInf[index][0] >= 0 and coordInf[index][0] < ImF.shape[0]:
                ImF[coordInf[index][0], coordInf[index][1], :] = [255,0,0]
                if coordInf[index][0]+1 >= 0 and coordInf[index][0]+1 < ImF.shape[0]:
                    ImF[coordInf[index][0]+1, coordInf[index][1], :] = [255,0,0]
                if coordInf[index][0]-1 >= 0 and coordInf[index][0]-1 < ImF.shape[0]:
                    ImF[coordInf[index][0]-1, coordInf[index][1], :] = [255,0,0]  
        
                
        #fascicles
        for a in range(len(splines_fasc)):
            newy = np.arange(-2*USimage.shape[1], 3*USimage.shape[1], 1)
            newx = np.int32(splines_fasc[a](newy))
            newy = newy + 2*USimage.shape[1]
            coord = np.vstack((newx, newy)).T
            for b in range(coord.shape[0]):
                if coord[b][0]>=0 and coord[b][0]>=intersecU[a][0]\
                and coord[b][0]<=intersecL[a][0] and coord[b][0]<ImF.shape[0]:
                    ImF[coord[b][0], coord[b][1], :] = [0,255,0]
                    if coord[b][1]-1>=0 and coord[b][1]-1<ImF.shape[1]:
                        ImF[coord[b][0], coord[b][1]-1, :] = [0,255,0]
                    if coord[b][1]+1>=0 and coord[b][1]+1<ImF.shape[1]:
                        ImF[coord[b][0], coord[b][1]+1, :] = [0,255,0]
    
    #                if coord[b][0]+1>=0 and coord[b][0]+1<ImF.shape[0]:
    #                    ImF[coord[b][0]+1, coord[b][1], :] = [0,255,0]
    #                if coord[b][0]-1>=0 and coord[b][0]-1<ImF.shape[0]:  
    #                    ImF[coord[b][0]-1, coord[b][1], :] = [0,255,0]
        
           
        #intersection points
        for i1 in range(len(splines_fasc)):
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
        cv2.imwrite(path_to_img[:-8]+'_final.jpg', ImF)
        cv2.imshow('Final image', ImF)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        
        '''
        #FEATURES
        import matplotlib.pyplot as plt
        fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharex = False, sharey = False, figsize =(15,15))
        
        #muscle thickness
        ax1.plot(absc, thickness, 'b+', markersize = 3)
        ax1.set_title('Muscle thickness')
        ax1.set_ylabel('thickness (mm)')
        ax1.set_xlabel('muscle axis (mm)')
        ax1.set_ybound(0,30)
        
        #pennation angle
        ax2.plot(PAu, PAd, 'g+', markersize = 5, label = 'inside image')
        ax2.plot(outPAu, outPAd, 'r+', markersize = 5, label = 'extrapolated')
        ax2.legend(loc = 'upper left', prop={'size': 6})
        ax2.set_title('Pennation angles')
        ax2.set_xlabel('Angle (degrees) / upper apo')
        ax2.set_ylabel('Angle (degrees) / deep apo')
        ax2.set_xbound(0,30)
        ax2.set_ybound(0,30)
    
        #fascicles length
        ax3.plot(fasc_in, [1] * len(fasc_in), 'g+', markersize = 5, label = 'inside image')
        ax3.plot(fasc_out, [1] * len(fasc_out), 'r+', markersize = 5, label = 'extrapolated')
        ax3.legend(loc = 'upper left', prop={'size': 6})
        ax3.set_title('Fascicles length')
        ax3.set_xlabel('length (mm)')
        
        plt.show()
        '''
        
        return archi_auto
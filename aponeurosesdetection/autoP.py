from calibration.calib import autoCalibration
from preprocessing.cropping import manualcropping
from preprocessing.preprocess import preprocessingApo
import apoLoc as apoL
import apoCont as apoC
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
#RGBimageP = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/jHamON_data/06_thomasmartine/fam_2/architecture/Martine_thomas_20180714_095546_image_bfp.jpg', -1)
#RGBimageP = cv2.imread('C:/Users/Lisa Paillard/Desktop/Pour TFE/jHamON_data/12_sufyan/fam_1/architecture/Sufyan_jamon_20180709_125927_image_bfp.jpg', -1)
#USimageP, pt_intersection = manualcropping(RGBimageP, 'C:/Users/Lisa Paillard/Desktop/Pour TFE/jHamON_data/06_thomasmartine/fam_2/architecture/Martine_thomas_20180714_095546_image_bfp.txt')
#USimageP, pt_intersection = manualcropping(RGBimageP, 'C:/Users/Lisa Paillard/Desktop/Pour TFE/jHamON_data/12_sufyan/fam_1/architecture/Sufyan_jamon_20180709_125927_image_bfp.txt')



def panoprocessing(path_to_image, path_to_txtfile):
    
    #opening image and validate processing
    RGBimageP = cv2.imread(path_to_image, -1)
    cv2.imshow('Image to process', RGBimageP)
    process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
        Close all windows after clicking yes or no.', default = 'yes', icon='question')
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    if process == True:
        
        #################################################
        
        #Calibrate the image
        calibX, calibY = autoCalibration(RGBimageP)
        if calibX > 2 * calibY or calibY > 2 * calibX:
            calibX = min(calibX, calibY)
            calibY = min(calibX, calibY)
            
        #################################################
        
        #Crop the image thanks to manual labelling
        #pt_intersection is the point where aponeuroses meet = insertion point
    
        USimageP, pt_intersection, l1,l2,c1,c2 = manualcropping(RGBimageP, path_to_txtfile)

        cv2.imshow('Cropped image', USimageP)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        
        
        #################################################
        #resize
        #update calibration factors and pt_intersection's coordinates
        initialsize = (USimageP.shape[1], USimageP.shape[0])
        PERCENTAGE = 160
        newWidth = int(initialsize[0]*PERCENTAGE/100)
        newHeight = int(initialsize[1]*PERCENTAGE/100)
        USimageP = cv2.resize(src = USimageP, dsize = (newWidth, newHeight), interpolation = cv2.INTER_CUBIC)
        cv2.imwrite(path_to_image[:-8]+'_cropped.jpg', USimageP)

        calibX = calibX / PERCENTAGE * 100
        calibY = calibY / PERCENTAGE * 100
        
        pt_intersection = (pt_intersection[0] * PERCENTAGE / 100, pt_intersection[1] * PERCENTAGE / 100)
        
        #################################################
        
        #Sampling to analyze band-by-band the image
        #NBANDS: number of bands used to sample the image
        #MAXBAND: limit below which both aponeuroses are analysed (beyong only upper apo)
        sampleSize = USimageP.shape[0]
        NBANDS = int(USimageP.shape[1]/sampleSize)
        if NBANDS%2 == 0:
            MAXBAND = int(NBANDS / 2) 
        else:
            MAXBAND = int(NBANDS / 2) + 1
                # sample1 = USimageP[:,:sampleSize]
                # sample2 = USimageP[:,sampleSize:2*sampleSize]
                # etc
        
        #################################################
        
        #Preprocessing
        USimageP_pp = preprocessingApo(I = USimageP, typeI = 'panoramic', mode = 'localmean', margin = 0, sizeContrast = 41)
        cv2.imwrite(path_to_image[:-8]+'_preprocessed.jpg', USimageP_pp)
        '''
        cv2.imshow('Pre-processed image',USimageP_pp)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        '''
        
        #################################################   
        #location of aponeuroses and linear approximation of aponeuroses
        
        contoursSup = []
        contoursInf = []
        
        for i in range(NBANDS):
            #
            #
            if i < MAXBAND: #look for both deep and superficial aponeuroses
                
                
                #find approximate locations and linear approximations
                paramSup, paramInf, locSup, locInf = apoL.twoApoLocation(USimageP_pp[:, i*sampleSize:(i+1)*sampleSize], angle1 = 80, angle2 = 101, thresh = None, calibV = calibX)
                
                if paramSup[0] != 'error':

                    #######################
                    #######################
                    ### upper apo processing
                    #sub-images of superficial aponeurosis
                    Sup_i = np.copy(USimageP[locSup[0]:locSup[1], i*sampleSize:(i+1)*sampleSize])
                    Sup_i_pp = np.copy(USimageP_pp[locSup[0]:locSup[1], i*sampleSize:(i+1)*sampleSize]) #upper aponeurosis in sample i
            
                    #Initiate contour with linear approximation
                    iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'quadrangle_param', param = [paramSup[0], paramSup[1]-locSup[0], 10])       
                    
                    #Calculate contour with active contour model
                    contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Upper aponeurosis contour found in ', nSup_i, ' steps')
        
                    #verify contour has been detected and ask for validation
                    if np.amin(contourSup_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                        iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'quadrangle_param', param = [paramSup[0], paramSup[1] - locSup[0], 40])
                        contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                        print('Upper aponeurosis contour found in ', nSup_i, ' steps')
        
                    if np.amin(contourSup_i)<=0: #if the contour exists, extract it
                        contourSup_image_i, contourSup_points_i = apoC.extractContour(contourSup_i, Sup_i, offSetX = locSup[0], offSetY = i * sampleSize)
                        visu = np.copy(USimageP[:, i*sampleSize:(i+1)*sampleSize])
                        visu[locSup[0]:locSup[1],:] = contourSup_image_i
                        #ask for manual validation of the contour
                        cv2.imshow('Sample i', visu)
                        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                        cv2.waitKey(0) & 0xFF
                        cv2.destroyAllWindows()
                        
                        if valid == False:
                            #try a second time with a new initial contour
                            points = np.array([[0, 0], [int(Sup_i_pp.shape[0]/2), 0], [int(Sup_i_pp.shape[0]/2), Sup_i_pp.shape[1]], [0, Sup_i_pp.shape[1]]])
                            iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'set_of_points', setPoints = points)
                            contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                            print('Upper aponeurosis contour found in ', nSup_i, ' steps')
                            if np.amin(contourSup_i) <= 0 :
                                contourSup_image_i, contourSup_points_i = apoC.extractContour(contourSup_i, Sup_i, offSetX = locSup[0], offSetY = i * sampleSize)
                                visu = np.copy(USimageP[:, i*sampleSize:(i+1)*sampleSize])
                                visu[locSup[0]:locSup[1],:] = contourSup_image_i
                                #ask for manual validation of the contour
                                cv2.imshow('Sample i', visu)
                                valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                                cv2.waitKey(0) & 0xFF
                                cv2.destroyAllWindows()                
                        
                        
                        if valid == True:
                            #add contour_i to list 'contoursSup'
                            for elem in contourSup_points_i:
                                contoursSup.append(elem)
                    
                    
                    ####################
                    ####################
                    ### deep apo
                    #sub-images of deep aponeurosis
                    Inf_i = np.copy(USimageP[locInf[0]:locInf[1], i*sampleSize:(i+1)*sampleSize])
                    Inf_i_pp = np.copy(USimageP_pp[locInf[0]:locInf[1], i*sampleSize:(i+1)*sampleSize]) #deep aponeurosis in sample i
                    
                    #Initiate contour with linear approximation
                    iniInf_i = apoC.initiateContour(Inf_i_pp, typeC = 'quadrangle_param', param = [paramInf[0], paramInf[1]-locInf[0], 10])
        
                    #Calculate contour with active contour model
                    contourInf_i, nInf_i = apoC.activeContour(Inf_i_pp, iniInf_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Deep aponeurosis contour found in ', nInf_i, ' steps')
                    
                    #Verify contour has been detected and ask for validation
                    if np.amin(contourInf_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                        iniInf_i = apoC.initiateContour(Inf_i_pp, typeC = 'quadrangle_param', param = [paramInf[0], paramInf[1] - locInf[0], 40])
                        contourInf_i, nInf_i = apoC.activeContour(Inf_i_pp, iniInf_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                        print('Deep aponeurosis contour found in ', nInf_i, ' steps')
        
                    if np.amin(contourInf_i)<=0: #if the contour exists, extract it
                        contourInf_image_i, contourInf_points_i = apoC.extractContour(contourInf_i, Inf_i, offSetX = locInf[0], offSetY = i * sampleSize)
                        visu = np.copy(USimageP[:, i*sampleSize:(i+1)*sampleSize])
                        visu[locInf[0]:locInf[1],:] = contourInf_image_i
                        #ask for manual validation of the contour
                        cv2.imshow('Sample i', visu)
                        valid = tkbox.askyesno('Need user validation', 'Do you validate the contours? After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                        cv2.waitKey(0) & 0xFF
                        cv2.destroyAllWindows()
                        
                        if valid == False:
                            #try a second time with a new initial contour
                            points = np.array([[int(Inf_i_pp.shape[0]/2), 0], [Inf_i_pp.shape[0], 0], [Inf_i_pp.shape[0],Inf_i_pp.shape[1]], [int(Inf_i_pp.shape[0]/2), Inf_i_pp.shape[1]]])
                            iniInf_i = apoC.initiateContour(Inf_i_pp, typeC = 'set_of_points', setPoints = points)
                            contourInf_i, nInf_i = apoC.activeContour(Inf_i_pp, iniInf_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                            print('Deep aponeurosis contour found in ', nInf_i, ' steps')
                            if np.amin(contourInf_i) <= 0 :
                                contourInf_image_i, contourInf_points_i = apoC.extractContour(contourInf_i, Inf_i, offSetX = locInf[0], offSetY = i * sampleSize)
                                visu = np.copy(USimageP[:, i*sampleSize:(i+1)*sampleSize])
                                visu[locInf[0]:locInf[1],:] = contourInf_image_i
                                #ask for manual validation of the contour
                                cv2.imshow('Sample i', visu)
                                valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                                cv2.waitKey(0) & 0xFF
                                cv2.destroyAllWindows()                
                        
                        if valid == True:
                            #add contour_i to contourInf
                            for elem in contourInf_points_i:
                                contoursInf.append(elem)
                
            
            ##################
            ##################
            #look for superficial aponeurosis only
            elif i > MAXBAND and abs(i*sampleSize-min((i+1)*sampleSize, int(pt_intersection[1])))>USimageP.shape[0]/2:
                
                offS = 0
                offS2 = USimageP.shape[0]
                Sup_i = np.copy(USimageP[offS:offS2, i*sampleSize:min((i+1)*sampleSize, int(pt_intersection[1]))])
                Sup_i_pp = np.copy(USimageP_pp[offS:offS2, i*sampleSize:min((i+1)*sampleSize, int(pt_intersection[1]))])
    
                # find approximate location of aponeurosis
                param, loc = apoL.oneApoLocation(Sup_i_pp, thresh = None, calibV = calibX, angle1 = int(50), angle2 = int(90))
                
                if param[0] != 'error':
                    #initiate contour
                    iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'quadrangle_param', param = [param[0], param[1] - offS, 10])
        
                    #calculate contour with active contour model
                    contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                    print('Upper aponeurosis contour found in ', nSup_i, ' steps')
                    
                    if np.amin(contourSup_i) > 0: #try a second time with a bigger initial contour if no contour has been found
                        iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'quadrangle_param', param = [param[0], param[1] - offS, 40])
                        contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                        print('Upper aponeurosis contour found in ', nSup_i, ' steps')
                          
                    if np.amin(contourSup_i) <= 0:
                        contourSup_image_i, contourSup_points_i = apoC.extractContour(contourSup_i, Sup_i, offSetX = offS, offSetY = i * sampleSize)
                        #ask for manual validation of the contour
                        cv2.imshow('Upper aponeurosis contour on sample i', contourSup_image_i)
                        valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                        cv2.waitKey(0) & 0xFF
                        cv2.destroyAllWindows()
                        
                        if valid == False:
                            #try a second time with a new initial contour
                            points = np.array([[0, 0], [int(Sup_i_pp.shape[0]/2), 0],[int(Sup_i_pp.shape[0]/2), Sup_i_pp.shape[1]], [0, Sup_i_pp.shape[1]]])
                            iniSup_i = apoC.initiateContour(Sup_i_pp, typeC = 'set_of_points', setPoints = points)
                            contourSup_i, nSup_i = apoC.activeContour(Sup_i_pp, iniSup_i, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
                            print('Upper aponeurosis contour found in ', nSup_i, ' steps')
                            if np.amin(contourSup_i) <= 0 :
                                contourSup_image_i, contourSup_points_i = apoC.extractContour(contourSup_i, Sup_i, offSetX = offS, offSetY = i * sampleSize)
                                #ask for manual validation of the contour
                                cv2.imshow('Upper aponeurosis contour on sample i', contourSup_image_i)
                                valid = tkbox.askyesno('Need user validation', 'Do you validate the contour ? If no, this section will be ignored in the interpolation process. After clicking yes or no, please close the image windows to continue.', default = 'yes', icon='question')
                                cv2.waitKey(0) & 0xFF
                                cv2.destroyAllWindows()
                            
                        if valid == True:
                            #add contour_i to list 'contoursSup'
                            for elem in contourSup_points_i:
                                contoursSup.append(elem)
        
        contoursSup.append(pt_intersection)
        contoursInf.append(pt_intersection)
        
        if len(contoursSup) <=1 or len(contoursInf) <= 1:
            archi_auto = dict()
            archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
            archi_auto['calfct_to_mm before resize'] = {'vertical axis': calibX * PERCENTAGE / 100, 'horizontal axis': calibY* PERCENTAGE / 100}
            archi_auto['aposup'] = {'coords':'error'}
            archi_auto['apoinf'] = {'coords':'error'}
            return archi_auto
            
        
        #Interpolation and extrapolation
        spline_Sup = apoC.approximateApo(p = contoursSup, apoType = 'upper', I = USimageP, typeapprox = 'polyfit', d = 3)
        spline_Inf = apoC.approximateApo(p = contoursInf, apoType = 'lower', I = USimageP, typeapprox = 'polyfit', d = 2)
        
            
        #Calculate coordinates of aponeuroses
        coordSup, spline_Sup = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_Sup, interval = [0, int(pt_intersection[1])])
        coordInf, spline_Inf = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_Inf, interval = [0, int(pt_intersection[1])])
        
        
       #Fascicles detection
        all_snippets = []
        all_snippets_line = []
        
        for i in range(NBANDS-1):
            
            minRow = np.amax(coordSup[i*sampleSize:(i+1)*sampleSize,0])
            maxRow = np.amin(coordInf[i*sampleSize:(i+1)*sampleSize,0])
            ROI = USimageP[minRow:maxRow, i*sampleSize:min((i+1)*sampleSize, pt_intersection[1]), :]
            if ROI.size>0:
                #Enhance tube-like structures with MVEF method - Frangi -
                #let's consider that fascicle diameter is between 0.3 mm and 0.5 mm,
                #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
                print('MVEF running')
                sca = np.arange(round(0.3/calibX), round(0.7/calibX), 0.5)
                MVEF_image = FaDe.MVEF_2D(255-ROI, sca, [0.5, 0])
                cv2.imwrite(path_to_image[:-8]+'__'+str(i)+'_mvef.jpg', MVEF_image)
                
                '''
                cv2.imshow('MVEF',MVEF_image)
                cv2.waitKey(0) & 0xFF
                cv2.destroyAllWindows()
                '''
                #
                #threshold
                threshMVEF_percent = 85
                threshMVEF = np.percentile(MVEF_image, threshMVEF_percent)
                MVEF_image2 = cv2.threshold(MVEF_image, threshMVEF, 255, cv2.THRESH_BINARY)[1]
                cv2.imwrite(path_to_image[:-8]+'__'+str(i)+'_mvef2.jpg', MVEF_image2)
                
                #locate muscle snippets and filter them
                snippets, snippets_line = FaDe.locateSnippets(MVEF_image2, calibX, calibY,\
                                                              minLength = 5, \
                                                              offSetX = minRow, offSetY = i*sampleSize, im = USimageP)
                if snippets == 'error':
                    #move coords back to original RGB image coordinate system
                    coordSup[:,0] = np.int64(coordSup[:,0]*100/PERCENTAGE+l1)
                    coordSup[:,1] = np.int64(coordSup[:,1]*100/PERCENTAGE+c1)
                    coordInf[:,0] = np.int64(coordInf[:,0]*100/PERCENTAGE+l1)
                    coordInf[:,1] = np.int64(coordInf[:,1]*100/PERCENTAGE+c1)
                    miniS = np.amin(coordSup[:,1])
                    maxiS = np.amax(coordSup[:,1])
                    miniI = np.amin(coordInf[:,1])
                    maxiI = np.amax(coordInf[:,1])
                    mini = max(miniS, miniI)
                    maxi = min(maxiS, maxiI)
                    coordS = np.zeros((maxi-mini+1,2))
                    coordI = np.zeros((maxi-mini+1,2))
                    for index1 in range(mini, maxi+1):
                        pixelS = [coordSup[i,0] for i in range(coordSup.shape[0]) if coordSup[i,1] == index1]
                        pixelI = [coordInf[i,0] for i in range(coordInf.shape[0]) if coordInf[i,1] == index1]
                        coordS[index1-mini,0] = np.mean(pixelS)
                        coordS[index1-mini,1] = index1
                        coordI[index1-mini,0] = np.mean(pixelI)
                        coordI[index1-mini,1] = index1
                    archi_auto = dict()
                    archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
                    archi_auto['calfct_to_mm before resize'] = {'vertical axis': calibX * PERCENTAGE/100, 'horizontal axis': calibY*PERCENTAGE/100}
                    archi_auto['aposup'] = {'coords':coordS}
                    archi_auto['apoinf'] = {'coords':coordI}
                    #muscle thickness measurement
                    abscissa, thickness, thickness_spline = MUFeaM.muscleThickness(points1 = coordI, points2 = coordS, start = mini, end = maxi, calibV = calibX*PERCENTAGE/100, calibH = calibY*PERCENTAGE/100)
                    archi_auto['MT'] = {'coords': thickness, 'columns interval': [mini, maxi]}
                    return archi_auto
                    
                    
                all_snippets = all_snippets + snippets
                all_snippets_line = all_snippets_line + snippets_line
            
        fasc, fasc2 = FaDe.combineSnippets(USimageP, all_snippets, all_snippets_line, min_nb_sn = 3, thresh_alignment = 5)
    
        #transform the snippets, which are in fact the contours of the fascicles,
        #into lines by taking the mean of the contour, that is the mean on each
        #column of the contour. When branches exist, the mean is not computed (the
        #region is ignored)
    
        averages = [] 
        for i in range(len(fasc)):
            averages.append(FaDe.contourAverage(fasc[i]))
            
        averages2 = []
        for i in range(len(fasc2)):
            averages2.append(FaDe.contourAverage(fasc2[i]))
            
        #interpolations to get fascicles' curve
        splines_fasc = FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages, d = 2)
        splines_fasc = splines_fasc + FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages2, d = 1)

        #intersections of fascicles with aponeuroses (in pixels)
        intersecL, intersecU, splines_fasc = MUFeaM.findIntersections(spl_inf = spline_Inf, spl_sup = spline_Sup,\
                                                                      listSpl = splines_fasc, start = 0, search_interval=[-c1*PERCENTAGE/100, pt_intersection[1]])

        #Location of fascicles: (in mm from aponeuroses intersection point)
        loc_fasc = MUFeaM.locateFasc(intersecL, pt_intersection, calibY)
        
        #Pennation angles calculation (in degree)
        PA_Sup = MUFeaM.pennationAngles(spline_Sup, splines_fasc, intersecU, calibX, calibY)
        PA_Inf = MUFeaM.pennationAngles(spline_Inf, splines_fasc, intersecL, calibX, calibY)
        
        #Fascicles length calculation (in mm)
        fasc_length = MUFeaM.fasciclesLength(splines_fasc, intersecU, intersecL, calibX, calibY)
    
        #Dict containing info per image
        #move coords back to original RGB image coordinate system
        coordS = np.zeros(coordSup.shape)
        coordI = np.zeros(coordInf.shape)
        coordS[:,0] = np.int64(coordSup[:,0]*100/PERCENTAGE+l1)
        coordS[:,1] = np.int64(coordSup[:,1]*100/PERCENTAGE+c1)
        coordI[:,0] = np.int64(coordInf[:,0]*100/PERCENTAGE+l1)
        coordI[:,1] = np.int64(coordInf[:,1]*100/PERCENTAGE+c1)        
        miniS = np.amin(coordS[:,1])
        maxiS = np.amax(coordS[:,1])
        miniI = np.amin(coordI[:,1])
        maxiI = np.amax(coordI[:,1])
        mini = int(max(miniS, miniI))
        maxi = int(min(maxiS, maxiI))
        coordS2 = np.zeros((maxi-mini+1,2))
        coordI2 = np.zeros((maxi-mini+1,2))
        for index1 in range(mini, maxi+1):
            pixelS = [coordS[i,0] for i in range(coordS.shape[0]) if coordS[i,1] == index1]
            pixelI = [coordI[i,0] for i in range(coordI.shape[0]) if coordI[i,1] == index1]
            coordS2[index1-mini,0] = np.mean(pixelS)
            coordS2[index1-mini,1] = index1
            coordI2[index1-mini,0] = np.mean(pixelI)
            coordI2[index1-mini,1] = index1       
        archi_auto = dict()
        archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
        archi_auto['calfct_to_mm before resize'] = {'vertical axis': calibX * PERCENTAGE / 100, 'horizontal axis': calibY* PERCENTAGE / 100}
        archi_auto['aposup'] = {'coords':coordS2}
        archi_auto['apoinf'] = {'coords':coordI2}
        #muscle thickness measurement
        abscissa, thickness, thickness_spline = MUFeaM.muscleThickness(points1 = coordI2, points2 = coordS2, start = mini, end = maxi, calibV = calibX*PERCENTAGE/100, calibH = calibY*PERCENTAGE/100)
        archi_auto['MT'] = {'coords': thickness, 'columns interval': [mini, maxi]}
        
        for index in range(len(splines_fasc)):
            archi_auto['fsc_' + str(index+1)] = {}
            typeIU = 'out of image'
            typeIL = 'out of image'
            typeFL = 'out of image'
            intU = intersecU[index]
            intL = intersecL[index]
            if intU[0] >= 0 and intU[0] < USimageP.shape[0] and\
                intU[1] >= 0 and intU[1] < USimageP.shape[1]:
                typeIU = 'in image'
            if intL[0] >= 0 and intL[0] < USimageP.shape[0] and\
                intL[1] >= 0 and intL[1] < USimageP.shape[1]:
                typeIL = 'in image'
            if typeIU == 'in image' and typeIL == 'in image':
                typeFL = 'in image'
                
            archi_auto['fsc_' + str(index+1)] = {'dist from insertion in mm': loc_fasc[index]}
            
            archi_auto['fsc_' + str(index+1)]['PAsup'] = {'in/out of the image': typeIU,
                                                          'value in degree' : PA_Sup[index],
                                                          'intersection with apo': [int(intU[0]*100/PERCENTAGE+l1), int(intU[1]*100/PERCENTAGE+c1)]}

            archi_auto['fsc_' + str(index+1)]['PAinf'] = {'in/out of the image': typeIL,
                                                          'value in degree' : PA_Inf[index],
                                                          'intersection with apo': [int(intL[0]*100/PERCENTAGE+l1), int(intL[1]*100/PERCENTAGE+c1)]}
            
            archi_auto['fsc_' + str(index+1)]['FL'] = {'in/out of the image':typeFL,
                                                       'length in mm': fasc_length[index]}
    
        #####
        #Visualization
        for index in range(coordSup.shape[0]):
            if coordSup[index][0] >= 0 and coordSup[index][0] < USimageP.shape[0]:
                USimageP[coordSup[index][0], coordSup[index][1], :] = [255, 0, 0]
                if coordSup[index][0] +1 >= 0 and coordSup[index][0]+1 < USimageP.shape[0]:
                    USimageP[coordSup[index][0]+1, coordSup[index][1],:] = [255, 0, 0]
                if coordSup[index][0]-1 >= 0 and coordSup[index][0]-1 < USimageP.shape[0]:
                    USimageP[coordSup[index][0]-1, coordSup[index][1],:] = [255, 0, 0]
            if coordInf[index][0] >= 0 and coordInf[index][0] < USimageP.shape[0]:
                USimageP[coordInf[index][0], coordInf[index][1], :] = [255, 0, 0]
                if coordInf[index][0]+1 >= 0 and coordInf[index][0]+1 < USimageP.shape[0]:
                    USimageP[coordInf[index][0]+1, coordInf[index][1], :] = [255, 0, 0]
                if coordInf[index][0]-1 >= 0 and coordInf[index][0]-1 < USimageP.shape[0]:
                    USimageP[coordInf[index][0]-1, coordInf[index][1], :] = [255, 0, 0]
      
        #snippets
        couleurs = [[255,0,0], [0,255,0], [0,0,255], [255,255,0],[255,0,255], [0,255,255],\
                    [100,200,0],[100,200,100], [50,200,0],[50,100,50], [255,100,0],\
                    [120,120,255], [255,80,80],[0,100,200], [0,100,80], [255,255,255],\
                    [120,120,120], [50,100,150],[100,50,150], [150,100,50], [50,150,100],
                    [100,150,50],[150,50,100],[12,75,255],[40,140,40]]
        for f in range(len(fasc)):
            for g in range(fasc[f].shape[0]):
                USimageP[fasc[f][g,0], fasc[f][g,1], :] = couleurs[f]
        for f in range(len(fasc2)):
            for g in range(fasc2[f].shape[0]):
                USimageP[fasc2[f][g,0], fasc2[f][g,1], :] = couleurs[-f]
                
                
        #fascicles
        for n1 in range(len(splines_fasc)):
            newy = np.arange(0, USimageP.shape[1], 1)
            newx = np.int32(splines_fasc[n1](newy))
            coord = np.vstack((newx, newy)).T
            for n2 in range(coord.shape[0]):
                if coord[n2][0]>=0 and coord[n2][0]<USimageP.shape[0]:
                    USimageP[coord[n2][0], coord[n2][1], :] = [0,255,0]
                    if coord[n2][1] - 1 >= 0 and coord[n2][1] - 1 < USimageP.shape[1]:
                        USimageP[coord[n2][0], coord[n2][1] - 1, :] = [0,255,0]
                    if coord[n2][1] + 1 >= 0 and coord[n2][1] + 1 < USimageP.shape[1]:
                        USimageP[coord[n2][0], coord[n2][1] + 1, :] = [0,255,0]

        #pointsintersection
        for i0 in range(len(intersecU)):
            intU = intersecU[i0]
            intUp = [int(intU[0]*100/PERCENTAGE+l1), int(intU[1]*100/PERCENTAGE+c1)]
            intL = intersecL[i0]
            intLow = [int(intL[0]*100/PERCENTAGE+l1), int(intL[1]*100/PERCENTAGE+c1)]
            for i1 in range(-5,6):
                if intUp[0]+i1>=0 and intUp[0]+i1<RGBimageP.shape[0] and\
                intUp[1]+i1>=0 and intUp[1]+i1< RGBimageP.shape[1]:
                    RGBimageP[intUp[0]+i1,intUp[1],:] = [0,0,255]
                    RGBimageP[intUp[0],intUp[1]+i1,:] = [0,0,255]
            for i2 in range(-7,8):
                if intLow[0]+i2>=0 and intLow[0]+i2<RGBimageP.shape[0] and\
                intLow[1]+i2>=0 and intLow[1]+i2< RGBimageP.shape[1]:
                    RGBimageP[intLow[0]+i2,intLow[1],:] = [0,255,255]
                    RGBimageP[intLow[0],intLow[1]+i2,:] = [0,255,255]
        
        cv2.imwrite(path_to_image[:-8]+'_final.jpg', USimageP)
        cv2.imshow('full image', USimageP)
        cv2.imshow('original image',RGBimageP)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()

        return archi_auto
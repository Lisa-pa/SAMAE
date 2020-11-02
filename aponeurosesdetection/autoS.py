"""Auto processing of simple/standard images"""

def simpleprocessing(path_to_img):
    """
    Function that realizes the (semi) automatic processing of simple/standard US images of muscles

    inputs
        path_to_image (string): path to the image + name image + extension

    outputs:
        a dictionary containing the analyzed architecture of the image.
    """

    from calibration.calib import autoCalibration
    from preprocessing.cropping import autocropping
    from preprocessing.preprocess import preprocessingApo
    import apoLoc as apoL
    import apoCont as apoC
    import MUFeaM as MUFeaM
    import FaDe as FaDe

    import cv2
    import numpy as np
    import tkinter.messagebox as tkbox

    # open the image and validate the start of the processing
    RGBimage = cv2.imread(path_to_img, -1)
    cv2.imshow('Image to process', RGBimage)
    process = tkbox.askyesno('Need user approval', 'Do you accept to process this image?\
        Close all windows after clicking yes or no.', default = 'yes', icon='question')
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    if process == True:
        

        # Calibrate the image
        # calibX, calibY are the calibration factors in the vertical and horizontal directions
        calibX, calibY = autoCalibration(RGBimage)
        #check that one calibration factor does not have an implausble value compared to the other
        if calibX > 2 * calibY or calibY > 2 * calibX:
            calibX = min(calibX, calibY)
            calibY = min(calibX, calibY)
                    
        #Crop the image to keep essential US data
        print('Try thresholds (10,15,12,25,2,6)')
        #   automatic try
        USimage, l1, l2, c1, c2 = autocropping(RGBimage, 10., 15., 12., 25., calibY, additionalCrop1 = 2., additionalCrop2 = 6.)
        if USimage.size<=0: #if the previous cropping did not work, try a second automatic try
            print('Try thresholds (6,15,6,25,0,0)')
            USimage, l1, l2, c1, c2 = autocropping(RGBimage, 6., 15., 6., 25., calibY, additionalCrop1 = 0, additionalCrop2 = 0)
        
        #--------ask for validation; if cropping does not suit the user, ask the user for thresholds in manual entry
        # you can try new thresholds maximum 5 times
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
                print('Cropping failed. Need new thresholds')
                entry1 = input("Enter integer betw. 0 and 255 for minimum threshold of columns' mean (recommended values: 10 or 6): ")
                entry2 = input("Enter integer betw. 0 and 255 for maximum threshold of columns' mean (recommended value: 15): ")
                entry3 = input("Enter integer betw. 0 and 255 for minimum threshold of raws' mean (recommended values: 12 or 6): ")
                entry4 = input("Enter integer betw. 0 and 255 for maximum threshold of raws' mean (recommended value: 25): ")
                entry5 = input("Optional additional cropping in mm at top of image:")
                entry6 = input("Optional additional cropping in mm at bottom of image:")
                THRESH1 = int(entry1)
                THRESH2 = int(entry2)
                THRESH3 = int(entry3)
                THRESH4 = int(entry4)
                CROP1 = int(entry5)
                CROP2 = int(entry6)
                print(f'You entered the following parameters: {THRESH1, THRESH2, THRESH3, THRESH4, CROP1, CROP2}')
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
        #               upper aponeurosis was spotted.
        # locInf:        lines in between which deep aponeurosis was spotted
        
        # preprocess image
        USimage_pp = preprocessingApo(I = USimage, typeI = 'simple', mode = 'localmean', margin = 0, sizeContrast = 41) #pre_processing
        #locate both aponeuroses in image + linear modeling of aponeuroses
        print('Looking for aponeuroses')
        paramSup, paramInf, locSup, locInf = apoL.twoApoLocation(USimage_pp, angle1 = 80, angle2 = 100, thresh = None, calibV = calibX)
        cv2.imwrite(path_to_img[:-8]+'_preprocessed.jpg', USimage_pp)        
        Bspl_sup = 0
        Bspl_inf = 0
        
        if paramSup[0] == 'error':
            archi_auto = dict()
            archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
            archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
            archi_auto['aposup'] = {'coords':'error'}
            archi_auto['apoinf'] = {'coords':'error'}
            print('The detection of aponeuroses failed')
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
        #initiate contour by quadrangle centered around linear modeling
        iniSupApo = apoC.initiateContour(SupApo_pp, typeC = 'quadrangle_param', param = [paramSup[0], paramSup[1] - locSup[0], 8])
        #evolve curve with active contour model
        contourSup, nSup = apoC.activeContour(SupApo_pp,iniSupApo,0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
        print('Upper aponeurosis contour found in ', nSup, ' steps')
        #check contour has been found and ask for MANUAL validation from the user
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
            visu = np.copy(USimage)
            visu[locSup[0]:locSup[1],:] = contourSup_image
            cv2.imshow('Upper aponeurosis contour', visu)
            valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
            cv2.waitKey(0) & 0xFF
            cv2.destroyAllWindows()
            if valid == False: #not validated by the user, try again with new initial contour
                iniSupApo = apoC.initiateContour(SupApo_pp, typeC = 'quadrangle_param', param = [0, int(SupApo.shape[0]/2), 40])
                contourSup, nSup = apoC.activeContour(SupApo_pp, iniSupApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
                print('Upper aponeurosis contour found in ', nSup, ' steps')
                if np.min(contourSup)<=0:
                    contourSup_image, contourSup_points = apoC.extractContour(contourSup, SupApo, offSetX = locSup[0], offSetY = 0)
                    visu = np.copy(USimage)
                    visu[locSup[0]:locSup[1],:] = contourSup_image
                    cv2.imshow('Upper aponeurosis contour', visu)
                    valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
                    cv2.waitKey(0) & 0xFF
                    cv2.destroyAllWindows()                
                    
            if valid == True:   #B-spline to approximate aponeurosis if contour suits
                                #replace paramSup coefficients by the spline
                type_approx_UA = 'spline'
                #polynomial approximation:
                paramSup = apoC.approximateApo(p = contourSup_points, apoType = 'upper', I = SupApo_pp, typeapprox = 'polyfit', d = 1)
                #real detected contour: (to estimate MT)
                Bspl_sup = apoC.approximateApo(p = contourSup_points, apoType = 'upper', I = SupApo_pp, typeapprox = 'Bspline', d = 1)
            elif valid == False: #use linear approximation from radon transform
                type_approx_UA = 'linear'
    
        #initiate contour of deep aponeurosis with quandrangle centered around linear modeling
        ini_InfApo = apoC.initiateContour(InfApo_pp, typeC = 'quadrangle_param', param = [paramInf[0], paramInf[1] - locInf[0], 8])
        # evolve contour with active contour model
        contourInf, nInf = apoC.activeContour(InfApo_pp, ini_InfApo, 0.5, 0.01, 0.02, 3.0, 1.0, 1.0, 65.025, 0.10)
        print('Lower aponeurosis contour found in ', nInf, ' steps')
        #check that a contour has been detected and ask the user to validate
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
            visu = np.copy(USimage)
            visu[locInf[0]:locInf[1],:] = contourInf_image
            cv2.imshow('Lower aponeurosis contour', visu)
            valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
            cv2.waitKey(0) & 0xFF
            cv2.destroyAllWindows()   
            if valid == False:  #not validated by the user, try again with new initial contour
                ini_InfApo = apoC.initiateContour(InfApo_pp, typeC = 'quadrangle_param', param = [0, int(InfApo.shape[0]/2), 40])
                contourInf, nInf = apoC.activeContour(InfApo_pp, ini_InfApo, 0.3,0.01,0.02,3.0, 1.0, 1.0, 65.025, 0.10)
                print('Lower aponeurosis contour found in ', nInf, ' steps')
                if np.min(contourInf)<=0:
                    contourInf_image, contourInf_Points = apoC.extractContour(contourInf, InfApo, offSetX = locInf[0], offSetY = 0)
                    visu = np.copy(USimage)
                    visu[locInf[0]:locInf[1],:] = contourInf_image
                    cv2.imshow('Lower aponeurosis contour', visu)
                    valid = tkbox.askyesno('Need user validation', 'Do you validate the contour? If no, linear approximation will be used in the rest of the process. After clicking yes or no, please close the image windows to continue.', default = 'no', icon='question')
                    cv2.waitKey(0) & 0xFF
                    cv2.destroyAllWindows()    

            if valid == True:   #B-spline to approximate aponeurosis if contour suits
                                #replace paramInf coefficients by the spline
                type_approx_LA = 'spline'
                #polynomial approximation:
                paramInf = apoC.approximateApo(p = contourInf_Points, apoType = 'lower', I = InfApo_pp, typeapprox = 'polyfit', d = 1)
                #real detected contour: (to estimate MT)
                Bspl_inf = apoC.approximateApo(p = contourInf_Points, apoType = 'lower', I = InfApo_pp, typeapprox = 'Bspline', d = 1)
            elif valid == False:#else linear approximation from Radon transform is used
                type_approx_LA = 'linear'
        
    
        #calculate coordinates of aponeuroses
        #if aponeusosis linear, param is transformed from a list of parameters into 
        # a spline during function pointsCoordinates
        #       approximation of aponeuroses:
        approx_sup, paramSup = MUFeaM.pointsCoordinates(typeA = type_approx_UA, param = paramSup, interval = [0, USimage.shape[1]])
        approx_inf, paramInf = MUFeaM.pointsCoordinates(typeA = type_approx_LA, param = paramInf, interval = [0, USimage.shape[1]])
        #       pixels really detected (= before approximation):
        if Bspl_sup != 0:
            coordSup, Bspl_sup = MUFeaM.pointsCoordinates(typeA = 'spline', param = Bspl_sup, interval = [0, RGBimage.shape[1]-1-c1])
        else:
            coordSup, Bspl_sup = MUFeaM.pointsCoordinates(typeA = 'spline', param = paramSup, interval = [0, RGBimage.shape[1]-1-c1])
        if Bspl_inf != 0:
            coordInf, Bspl_inf = MUFeaM.pointsCoordinates(typeA = 'spline', param = Bspl_inf, interval = [0, RGBimage.shape[1]-1-c1])
        else:
            coordInf, Bspl_inf = MUFeaM.pointsCoordinates(typeA = 'spline', param = paramInf, interval = [0, RGBimage.shape[1]-1-c1])
        
        
        #Muscle thickness calculation
        absc, thickness, spline_thickness = MUFeaM.muscleThickness(points1 = coordSup, points2 = coordInf,\
                                                                   start = 0, end = RGBimage.shape[1]-1-c1, calibV = calibX, calibH = calibY)
            
        
        #Start fascicles detection
        print('Looking for muscle fascicles')
        #ROI : region of interest in between aponeuroses
        # contains only muscle fascicles
        crop1 = np.amax(approx_sup[:,0]) + 30
        crop2 = np.amin(approx_inf[:,0]) - 30
        ROI = np.copy(USimage[crop1:crop2, :, :])
        cv2.imwrite(path_to_img[:-8]+'_ROI.jpg', ROI)
                
        #Enhance tube-like structures with MVEF method - Frangi - 
        #let's consider that fascicle diameter is between 0.3 mm and 0.5 mm,
        #the following list is the equivalent interval in pixels, with a step of 0.5 pixel
        sca = np.arange(round(0.3/calibX), round(0.5/calibX), 0.5)
        MVEF_image = FaDe.MVEF_2D(255-ROI, sca, [0.5, 0])
        cv2.imwrite(path_to_img[:-8]+'_MVEF.jpg', MVEF_image)

        #threshold to binarize filtered image
        threshMVEF_percent = 85
        threshMVEF = np.percentile(MVEF_image, threshMVEF_percent)
        MVEF_image2 = cv2.threshold(MVEF_image, threshMVEF, 255, cv2.THRESH_BINARY)[1]
        cv2.imwrite(path_to_img[:-8]+'_MVEF2.jpg', MVEF_image2)
    
        #locate snippets and filter them
        snippets, snip_lines = FaDe.locateSnippets(MVEF_image2, calibX, calibY,\
                                                   minLength = 4,\
                                                   offSetX = crop1)
        
        if snippets == 'error':
            #move coords to original image coordinate system
            approx_sup[:,0] = approx_sup[:,0]+l1
            approx_sup[:,1] = approx_sup[:,1]+c1
            approx_inf[:,0] = approx_inf[:,0]+l1
            approx_inf[:,1] = approx_inf[:,1]+c1
            coordSup[:,0] = coordSup[:,0]+l1
            coordSup[:,1] = coordSup[:,1]+c1
            coordInf[:,0] = coordInf[:,0]+l1
            coordInf[:,1] = coordInf[:,1]+c1
            archi_auto = dict()
            archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
            archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
            archi_auto['aposup'] = {'coords':coordSup, 'approx':approx_sup}
            archi_auto['apoinf'] = {'coords':coordInf, 'approx':approx_inf}
            archi_auto['MT'] = {'coords': thickness, 'columns interval': [0+c1, RGBimage.shape[1]-1]}
            print('The detection of muscle fascicles failed.')
            return archi_auto 
        
        #look for aligned snippets, that would be part of the same fascicle
        #fasc: fascicles made of 2 snippets and more
        #fasc2: fascicles made of less than 2 snippets
        fasc, fasc2 = FaDe.combineSnippets(USimage, snippets, snip_lines, min_nb_sn = 2, thresh_alignment = 20)
    
        #transform the snippets, which are in fact the contours of the fascicles,
        #into lines
        averages = [] 
        for i in range(len(fasc)):
            averages.append(FaDe.contourAverage(fasc[i]))

        averages2 = [] 
        for i in range(len(fasc2)):
            averages2.append(FaDe.contourAverage(fasc2[i]))        
        
        
        #interpolations to get fascicles' curve
        #all with degree = 1
        splines_fasc = FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages, d = 1)
        splines_fasc = splines_fasc + FaDe.approximateFasc(typeapprox = 'polyfit', listF = averages2, d = 1)
        
        #intersections of fascicles with aponeuroses
        #determination of the dominant orientation of fascicles (positive or negative slope)
        #to eliminate fascicles that would not follow the dominant tendancy
        counter_posit = 0
        counter_negat = 0
        for fasci in splines_fasc:
            if fasci((c1+c2)/2 + 50)-fasci((c1+c2)/2 - 50) > 0:
                counter_posit =counter_posit+1
            elif fasci((c1+c2)/2 + 50)-fasci((c1+c2)/2 - 50) < 0:
                counter_negat =counter_negat+1
        if counter_negat >counter_posit:
            sig = -1
        if counter_negat< counter_posit:
            sig = 1
        intersecL, intersecU, splines_fasc = MUFeaM.findIntersections(spl_inf = paramInf, spl_sup = paramSup,\
                                                                      listSpl = splines_fasc, signOfSlope = sig, start = 0, search_interval = [-200/calibY, 200/calibY])
        
        print('Computing muscle architecture parameters.')
        #Pennation angles calculation (in degrees)
        PASup = MUFeaM.pennationAngles(paramSup, splines_fasc, intersecU, calibX, calibY)
        PA_Inf = MUFeaM.pennationAngles(paramInf, splines_fasc, intersecL, calibX, calibY)
        
        #Fascicles length calculation (in millimeters)
        fasc_length = MUFeaM.fasciclesLength(splines_fasc, intersecU, intersecL, calibX, calibY)
        
        #fascicle location in the original RGB image
        loc_fasc = MUFeaM.locateFasc(intersecL, [-l1,-c1], calibY)
        
        #Dict containing info per image
        #move coords to original image coordinate system

        approxS = np.zeros(approx_sup.shape)
        approxI = np.zeros(approx_inf.shape)
        approxS[:,0] = approx_sup[:,0]+l1
        approxS[:,1] = approx_sup[:,1]+c1
        approxI[:,0] = approx_inf[:,0]+l1
        approxI[:,1] = approx_inf[:,1]+c1
        coordSup[:,0] = coordSup[:,0]+l1
        coordSup[:,1] = coordSup[:,1]+c1
        coordInf[:,0] = coordInf[:,0]+l1
        coordInf[:,1] = coordInf[:,1]+c1        
        archi_auto = dict()
        archi_auto['crop'] = {'lines': [l1,l2], 'columns': [c1,c2]}
        archi_auto['calfct_to_mm'] = {'vertical axis': calibX, 'horizontal axis': calibY}
        archi_auto['aposup'] = {'coords':coordSup, 'approx':approxS}
        archi_auto['apoinf'] = {'coords':coordInf, 'approx':approxI}
        archi_auto['MT'] = {'coords': thickness, 'columns interval': [0+c1, RGBimage.shape[1]-1]}
        
        for index in range(len(splines_fasc)):
            archi_auto['fsc_' + str(index+1)] = {}
            typeIU = 'out of image'
            typeIL = 'out of image'
            typeFL = 'out of image'
            intU = intersecU[index]
            intL = intersecL[index]
            if intU[0] >= 0 and intU[0] < USimage.shape[0] and\
                intU[1] >= 0 and intU[1] < USimage.shape[1]:
                typeIU = 'in image'
            if intL[0] >= 0 and intL[0] < USimage.shape[0] and\
                intL[1] >= 0 and intL[1] < USimage.shape[1]:
                typeIL = 'in image'
            if typeIU == 'in image' and typeIL == 'in image':
                typeFL = 'in image' 
                
            archi_auto['fsc_' + str(index+1)] = {'dist from (0,0) of RGB image, in mm': loc_fasc[index]}
            archi_auto['fsc_' + str(index+1)]['PAsup'] = {'in/out of the image': typeIU,
                                                          'value in degree' : PASup[index],
                                                          'intersection with apo': [intU[0]+l1, intU[1]+c1]}

            archi_auto['fsc_' + str(index+1)]['PAinf'] = {'in/out of the image': typeIL,
                                                          'value in degree' : PA_Inf[index],
                                                          'intersection with apo': [intL[0]+l1, intL[1]+c1]}
            
            archi_auto['fsc_' + str(index+1)]['FL'] = {'in/out of the image':typeFL,
                                                       'length in mm': fasc_length[index]}   
    
    
    
    
    
    
    
    
    
    
    
        ################
        #Visualization of modeled aponeuroses, detected snippets and modeled fascicles
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
        approx_sup = np.vstack((newx_u, newy)).T
        approx_inf = np.vstack((newx_l, newy)).T
        
        for index in range(approx_sup.shape[0]):
            if approx_sup[index][0] >= 0 and approx_sup[index][0] < ImF.shape[0]:
                ImF[approx_sup[index][0], approx_sup[index][1],:] = [255,0,0]
                if approx_sup[index][0] +1 >= 0 and approx_sup[index][0]+1 < ImF.shape[0]:
                    ImF[approx_sup[index][0]+1, approx_sup[index][1],:] = [255,0,0]
                if approx_sup[index][0]-1 >= 0 and approx_sup[index][0]-1 < ImF.shape[0]:
                    ImF[approx_sup[index][0]-1, approx_sup[index][1],:] = [255,0,0]
            if approx_inf[index][0] >= 0 and approx_inf[index][0] < ImF.shape[0]:
                ImF[approx_inf[index][0], approx_inf[index][1], :] = [255,0,0]
                if approx_inf[index][0]+1 >= 0 and approx_inf[index][0]+1 < ImF.shape[0]:
                    ImF[approx_inf[index][0]+1, approx_inf[index][1], :] = [255,0,0]
                if approx_inf[index][0]-1 >= 0 and approx_inf[index][0]-1 < ImF.shape[0]:
                    ImF[approx_inf[index][0]-1, approx_inf[index][1], :] = [255,0,0]  
        
                
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

        ImF = cv2.resize(src = ImF, dsize = (int(ImF.shape[1]*0.4), int(ImF.shape[0]*0.4)), interpolation = cv2.INTER_CUBIC)
        cv2.imwrite(path_to_img[:-8]+'_final.jpg', ImF)
        cv2.imshow('Final image', ImF)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
        
        return archi_auto
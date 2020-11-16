import matplotlib.pyplot as plt
from SAMAE.plot import rand_jitter

def plotPairingSimple(path_to_dict, name_dict, colors, participants):
    """Function that aims at comparing the way fascicles should be matched
    in simple images (in particular, maximum thresholds between match auto and manual fascicles)
    path_to_dict: path to folder containing the dictionnary
    name_dict: name of the dictionnary within the folder
    colors: list of color tuples; length should be at least equal to
    the number of analysed participants.
    participant: list of analysed participants' names/ID.     
    """
    import numpy as np
    import fnmatch
    from dictmanager import load_obj
    
    #d = list of distances from reference point
    #fl = list of FL
    #PAs = list of PA sup
    #PAi = list of PA inf
    #mt = list of MT
    #_s : simple images
    #_m : manual
    #_a : automated

    '************************************************************************'
    '*****************************INITIALIZATION*****************************'
    
    spacial_thresh = [20,15,10,5,3,2,1,0] #millimeters    
    
    d_s_m = [[] for par in range(len(participants))]
    fl_s_m = [[] for par in range(len(participants))]
    PAs_s_m = [[] for par in range(len(participants))]
    PAi_s_m = [[] for par in range(len(participants))]
    mt_s_m = [[] for par in range(len(participants))]
    d_s_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    fl_s_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAs_s_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAi_s_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]

    d_s_a = [[] for par in range(len(participants))]
    fl_s_a = [[] for par in range(len(participants))]
    PAs_s_a = [[] for par in range(len(participants))]
    PAi_s_a = [[] for par in range(len(participants))]
    mt_s_a = [[] for par in range(len(participants))]
    d_s_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    fl_s_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAs_s_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAi_s_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
 
    diff_calfct_s = []



    #plots to determine how to pair auto and manual fascicles
    manuPairing_d_s_m=[[] for par in range(len(participants))]
    manuPairing_d_s_a=[[] for par in range(len(participants))]
    manuPairing_fl_s_m=[[] for par in range(len(participants))]
    manuPairing_fl_s_a=[[] for par in range(len(participants))]
    manuPairing_pas_s_m=[[] for par in range(len(participants))]
    manuPairing_pas_s_a=[[] for par in range(len(participants))]
    manuPairing_pai_s_m=[[] for par in range(len(participants))]
    manuPairing_pai_s_a=[[] for par in range(len(participants))]
    d_s_m_rel=[]
    d_s_a_rel=[]
    fl_s_m_rel=[]
    fl_s_a_rel=[]
    pas_s_m_rel=[]
    pas_s_a_rel=[]
    pai_s_m_rel=[]
    pai_s_a_rel=[]

    figSD, axSD = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figSD.suptitle('Distance to the distal insertion point (mm) depending on the pairing, in simple images', fontsize=40)
    axSD.set_ylabel('Distance to distal insertion (in mm)', fontsize= 35)
    axSD.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axSD.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axSD.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axSD.grid(True)
    axSD.tick_params(axis='x', labelsize=28, colors='k')
    axSD.tick_params(axis='y', labelsize=25, colors='k')

    figSfl, axSfl = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figSfl.suptitle('FL (mm) depending on the pairing, in simple images', fontsize=40)
    axSfl.set_ylabel('FL (in mm)', fontsize= 35)
    axSfl.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axSfl.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axSfl.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axSfl.grid(True)
    axSfl.tick_params(axis='x', labelsize=28, colors='k')
    axSfl.tick_params(axis='y', labelsize=25, colors='k')

    figSpas, axSpas = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figSpas.suptitle('PA superior (degree) depending on the pairing, in simple images', fontsize=40)
    axSpas.set_ylabel('PA superior (degree)', fontsize= 35)
    axSpas.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axSpas.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axSpas.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axSpas.grid(True)
    axSpas.tick_params(axis='x', labelsize=28, colors='k')
    axSpas.tick_params(axis='y', labelsize=25, colors='k')

    figSpai, axSpai = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figSpai.suptitle('PA inferior (degree) depending on the pairing, in simple images', fontsize=40)
    axSpai.set_ylabel('PA inferior (degree)', fontsize= 35)
    axSpai.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axSpai.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axSpai.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axSpai.grid(True)
    axSpai.tick_params(axis='x', labelsize=28, colors='k')
    axSpai.tick_params(axis='y', labelsize=25, colors='k')

    
    figSrel, axSrel = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figSrel.suptitle('Influence of the distance from reference point on fascicle length estimation (simple images)', fontsize=40)
    axSrel.set_xlabel('Distance from reference point (in mm)', fontsize= 35)
    axSrel.set_ylabel('Manual and auto fascicles length (in mm)', fontsize= 35)
    axSrel.grid(True)
    axSrel.tick_params(axis='x', labelsize=28, colors='k')
    axSrel.tick_params(axis='y', labelsize=25, colors='k')


    '************************************************************************'
    '*****************************DATA RETRIEVAL*****************************'    
    dictio = load_obj(name_dict, path_to_dict)
    l2 = ['fasc*', 'fsc_*']
    
    for par in range(len(participants)):
        
        participant = participants[par]
        fam_folders = [str(d) for d in dictio[participant].keys()]

        s_manuFasc = []
        s_autoFasc = []
        
        for fam in fam_folders:
            
            ###################################################################
            # simple images
            dictioS = dictio[participant][fam]['BF']['simple']
            images = [str(im) for im in dictioS.keys()]
            for i in images:
                
                ###############################################################
                # SIMPLE - manual
                dictioM = dictioS[i]['architecture manual']
                fascicles = [str(fa) for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    idf = fam + '/' + i + '/' + f
                    if len(dictioF.keys())>1:
                        s_manuFasc.append((idf, dictioF['dist from (0,0) of RGB image, in mm']))
                        d_s_m[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                        fl_s_m[par].append(dictioF['FL']['length in mm'])
                        PAs_s_m[par].append(dictioF['PAsup']['value in degree'])
                        PAi_s_m[par].append(dictioF['PAinf']['value in degree'])

                ###############################################################
                # SIMPLE - automatic
                if ('architecture auto' in dictioS[i]):
                    dictioA = dictioS[i]['architecture auto']
                    midRow = np.mean(dictioA['crop']['lines'])
                    midCol = np.mean(dictioA['crop']['columns'])
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        for f in fascicles:
                            dictioF = dictioA[f]
                            idf = fam + '/' + i + '/' + f
                            if len(dictioF.keys())>1:
                                #keep the fascicles that are in the lower half of the image,
                                #to compare with manual data - often taken in that region
                                PAi = dictioF['PAinf']['intersection with apo']
                                PAs = dictioF['PAsup']['intersection with apo']
                                fasc_row = (PAs[0]-PAi[0])/(PAs[1]-PAi[1])*(midCol-PAs[1])+PAs[0]
                                
                                if fasc_row <= midRow:
                                    s_autoFasc.append((idf, dictioF['dist from (0,0) of RGB image, in mm']))
                                    d_s_a[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                                    fl_s_a[par].append(dictioF['FL']['length in mm'])
                                    PAs_s_a[par].append(dictioF['PAsup']['value in degree'])
                                    PAi_s_a[par].append(dictioF['PAinf']['value in degree'])
                                    
                                
                        if ('MT for labelled points' in dictioM['MT']):
                            for ind0 in range(len(dictioM['MT']['MT for labelled points'])):
                                elem = dictioM['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_s_m[par].append(elem) #MT in mm
                            
                            for ind0 in range(len(dictioA['MT']['MT for labelled points'])):
                                elem = dictioA['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_s_a[par].append(elem)
                                            
                    
                    #calibration factors difference
                    diff_calfct_s.append(abs(dictioM['calfct_to_mm']-\
                                             dictioA['calfct_to_mm']['horizontal axis']))
            

        '************************************************************************'
        '***************SEARCH FOR MATCHING AUTO & MANUAL FASCICLES**************'
        listePair_manuF_s = []

        for n in range(len(s_manuFasc)):
            mf = s_manuFasc[n]
            subtr = [(tup,abs(tup[1]- mf[1])) for tup in s_autoFasc]
            subtr.sort(key=lambda x:x[1])
            closest = subtr[0]
            listePair_manuF_s.append((mf[0], closest[0][0], closest[1])) #tuple = ( ID manu fasc, ID auto fasc, distance entre les deux)
        listePair_manuF_s.sort(key=lambda x:x[1])
        uniqueMatching = []
        counterL = 0
        while counterL < len(listePair_manuF_s):
            currentAutoFasc = listePair_manuF_s[counterL][1]
            correspondingAutoFasc = [(listePair_manuF_s[counterL][0], listePair_manuF_s[counterL][2])]
            rank = counterL + 1
            while rank<len(listePair_manuF_s) and listePair_manuF_s[rank][1] == currentAutoFasc:
                correspondingAutoFasc.append((listePair_manuF_s[rank][0],listePair_manuF_s[rank][2]))
                rank = rank + 1
            correspondingAutoFasc.sort(key=lambda x:x[1])
            uniqueMatching.append((correspondingAutoFasc[0][0], currentAutoFasc, correspondingAutoFasc[0][1]))
            counterL = rank

        
        #manual fascicles pairing with an auto fascicle (manual fascicles can be used for several auto fasc)  
        for el in listePair_manuF_s:
            pathA0 = el[1].split('/')
            pathM0 = el[0].split('/')
            manuPairing_d_s_m[par].append(dictio[participant][pathM0[0]]['BF']['simple'][pathM0[1]]['architecture manual'][pathM0[2]]['dist from (0,0) of RGB image, in mm'])
            manuPairing_d_s_a[par].append(dictio[participant][pathA0[0]]['BF']['simple'][pathA0[1]]['architecture auto'][pathA0[2]]['dist from (0,0) of RGB image, in mm'])
            manuPairing_fl_s_m[par].append(dictio[participant][pathM0[0]]['BF']['simple'][pathM0[1]]['architecture manual'][pathM0[2]]['FL']['length in mm'])
            manuPairing_fl_s_a[par].append(dictio[participant][pathA0[0]]['BF']['simple'][pathA0[1]]['architecture auto'][pathA0[2]]['FL']['length in mm'])
            manuPairing_pas_s_m[par].append(dictio[participant][pathM0[0]]['BF']['simple'][pathM0[1]]['architecture manual'][pathM0[2]]['PAsup']['value in degree'])
            manuPairing_pas_s_a[par].append(dictio[participant][pathA0[0]]['BF']['simple'][pathA0[1]]['architecture auto'][pathA0[2]]['PAsup']['value in degree'])
            manuPairing_pai_s_m[par].append(dictio[participant][pathM0[0]]['BF']['simple'][pathM0[1]]['architecture manual'][pathM0[2]]['PAinf']['value in degree'])
            manuPairing_pai_s_a[par].append(dictio[participant][pathA0[0]]['BF']['simple'][pathA0[1]]['architecture auto'][pathA0[2]]['PAinf']['value in degree'])
        
        #reciprocal matching without any threshold on the distance between auto and manual fascicles
        for element in uniqueMatching:
            pathA1 = element[1].split('/')
            pathM1 = element[0].split('/')
            d_s_m_rel.append(dictio[participant][pathM1[0]]['BF']['simple'][pathM1[1]]['architecture manual'][pathM1[2]]['dist from (0,0) of RGB image, in mm'])
            d_s_a_rel.append(dictio[participant][pathA1[0]]['BF']['simple'][pathA1[1]]['architecture auto'][pathA1[2]]['dist from (0,0) of RGB image, in mm'])
            fl_s_m_rel.append(dictio[participant][pathM1[0]]['BF']['simple'][pathM1[1]]['architecture manual'][pathM1[2]]['FL']['length in mm'])
            fl_s_a_rel.append(dictio[participant][pathA1[0]]['BF']['simple'][pathA1[1]]['architecture auto'][pathA1[2]]['FL']['length in mm'])
            pas_s_m_rel.append(dictio[participant][pathM1[0]]['BF']['simple'][pathM1[1]]['architecture manual'][pathM1[2]]['PAsup']['value in degree'])
            pas_s_a_rel.append(dictio[participant][pathA1[0]]['BF']['simple'][pathA1[1]]['architecture auto'][pathA1[2]]['PAsup']['value in degree'])
            pai_s_m_rel.append(dictio[participant][pathM1[0]]['BF']['simple'][pathM1[1]]['architecture manual'][pathM1[2]]['PAinf']['value in degree'])
            pai_s_a_rel.append(dictio[participant][pathA1[0]]['BF']['simple'][pathA1[1]]['architecture auto'][pathA1[2]]['PAinf']['value in degree'])
        
        
        #reciprocal pairing with threshold th
        for t in range(len(spacial_thresh)):
            th = spacial_thresh[t]
            for element in uniqueMatching:
                pathA = element[1].split('/')
                pathM = element[0].split('/')
                
                if element[2] < th :
    
                    d_s_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['dist from (0,0) of RGB image, in mm'])
                    fl_s_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['FL']['length in mm'])
                    PAs_s_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['PAsup']['value in degree'])
                    PAi_s_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['PAinf']['value in degree'])
                    
                    d_s_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['dist from (0,0) of RGB image, in mm'])
                    fl_s_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['FL']['length in mm'])
                    PAs_s_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['PAsup']['value in degree'])
                    PAi_s_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['PAinf']['value in degree'])

            #reciprocal pairing with threshold th
            axSD.plot(rand_jitter([2*t-0.4]*len(d_s_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_s_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axSD.plot(rand_jitter([2*t+0.4]*len(d_s_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_s_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axSfl.plot(rand_jitter([2*t-0.4]*len(fl_s_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_s_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axSfl.plot(rand_jitter([2*t+0.4]*len(fl_s_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_s_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axSpas.plot(rand_jitter([2*t-0.4]*len(PAs_s_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAs_s_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axSpas.plot(rand_jitter([2*t+0.4]*len(PAs_s_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAs_s_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axSpai.plot(rand_jitter([2*t-0.4]*len(PAi_s_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAi_s_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axSpai.plot(rand_jitter([2*t+0.4]*len(PAi_s_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAi_s_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

        #all manu fascicles were paired with the closest auto fasc -> one auto fasc could be used several time
        axSD.plot(rand_jitter([-4-0.4]*len(manuPairing_d_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_d_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSD.plot(rand_jitter([-4+0.4]*len(manuPairing_d_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_d_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSfl.plot(rand_jitter([-4-0.4]*len(manuPairing_fl_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_fl_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSfl.plot(rand_jitter([-4+0.4]*len(manuPairing_fl_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_fl_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSpas.plot(rand_jitter([-4-0.4]*len(manuPairing_pas_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_pas_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSpas.plot(rand_jitter([-4+0.4]*len(manuPairing_pas_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_pas_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSpai.plot(rand_jitter([-4-0.4]*len(manuPairing_pai_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_pai_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSpai.plot(rand_jitter([-4+0.4]*len(manuPairing_pai_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_pai_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

        #all fascicles, not even paired
        axSD.plot(rand_jitter([-6-0.4]*len(d_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSD.plot(rand_jitter([-6+0.4]*len(d_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSfl.plot(rand_jitter([-6-0.4]*len(fl_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSfl.plot(rand_jitter([-6+0.4]*len(fl_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSpas.plot(rand_jitter([-6-0.4]*len(PAs_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAs_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSpas.plot(rand_jitter([-6+0.4]*len(PAs_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAs_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axSpai.plot(rand_jitter([-6-0.4]*len(PAi_s_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAi_s_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axSpai.plot(rand_jitter([-6+0.4]*len(PAi_s_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAi_s_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

    for t2 in range(len(spacial_thresh)):
        
        mean_d_s_a_th =np.mean([item for sublist in d_s_a_filtered[t2] for item in sublist])
        mean_d_s_m_th = np.mean([item for sublist in d_s_m_filtered[t2] for item in sublist])
        mean_fl_s_a_th =np.mean([item for sublist in fl_s_a_filtered[t2] for item in sublist])
        mean_fl_s_m_th = np.mean([item for sublist in fl_s_m_filtered[t2] for item in sublist])
        mean_pas_s_m_th=np.mean([item for sublist in PAs_s_m_filtered[t2] for item in sublist])
        mean_pas_s_a_th=np.mean([item for sublist in PAs_s_a_filtered[t2] for item in sublist])
        mean_pai_s_m_th=np.mean([item for sublist in PAi_s_m_filtered[t2] for item in sublist])
        mean_pai_s_a_th=np.mean([item for sublist in PAi_s_a_filtered[t2] for item in sublist])
        axSD.plot(2*t2-0.4, mean_d_s_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSD.plot(2*t2+0.4, mean_d_s_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSfl.plot(2*t2-0.4, mean_fl_s_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSfl.plot(2*t2+0.4, mean_fl_s_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSpas.plot(2*t2-0.4, mean_pas_s_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSpas.plot(2*t2+0.4, mean_pas_s_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSpai.plot(2*t2-0.4, mean_pai_s_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axSpai.plot(2*t2+0.4, mean_pai_s_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
    
    mean_d_s_m = np.mean([item for sublist in d_s_m for item in sublist])
    mean_d_s_a = np.mean([item for sublist in d_s_a for item in sublist])
    mean_fl_s_m = np.mean([item for sublist in fl_s_m for item in sublist])
    mean_fl_s_a = np.mean([item for sublist in fl_s_a for item in sublist])
    mean_pas_s_m = np.mean([item for sublist in PAs_s_m for item in sublist])
    mean_pas_s_a =np.mean([item for sublist in PAs_s_a for item in sublist])
    mean_pai_s_m=np.mean([item for sublist in PAi_s_m for item in sublist])
    mean_pai_s_a=np.mean([item for sublist in PAi_s_a for item in sublist])
    meanManuPair_d_s_m = np.mean([item for sublist in manuPairing_d_s_m for item in sublist])
    meanManuPair_d_s_a = np.mean([item for sublist in manuPairing_d_s_a for item in sublist])
    meanManuPair_fl_s_m = np.mean([item for sublist in manuPairing_fl_s_m for item in sublist])
    meanManuPair_fl_s_a = np.mean([item for sublist in manuPairing_fl_s_a for item in sublist])
    meanManuPair_pas_s_m = np.mean([item for sublist in manuPairing_pas_s_m for item in sublist])
    meanManuPair_pas_s_a = np.mean([item for sublist in manuPairing_pas_s_a for item in sublist])
    meanManuPair_pai_s_m = np.mean([item for sublist in manuPairing_pai_s_m for item in sublist])
    meanManuPair_pai_s_a = np.mean([item for sublist in manuPairing_pai_s_a for item in sublist])
    meanUniqueMatch_d_s_m = np.mean(d_s_m_rel)
    meanUniqueMatch_d_s_a = np.mean(d_s_a_rel)
    meanUniqueMatch_fl_s_m = np.mean(fl_s_m_rel)
    meanUniqueMatch_fl_s_a = np.mean(fl_s_a_rel)
    meanUniqueMatch_pas_s_m = np.mean(pas_s_m_rel)
    meanUniqueMatch_pas_s_a = np.mean(pas_s_a_rel)
    meanUniqueMatch_pai_s_m = np.mean(pai_s_m_rel)
    meanUniqueMatch_pai_s_a = np.mean(pai_s_a_rel)

    axSD.plot(-6.4, mean_d_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSD.plot(-5.6, mean_d_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSD.plot(-4.4, meanManuPair_d_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSD.plot(-3.6, meanManuPair_d_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot(-6.4, mean_fl_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot(-5.6, mean_fl_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot(-4.4, meanManuPair_fl_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot(-3.6, meanManuPair_fl_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot(-6.4, mean_pas_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot(-5.6, mean_pas_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot(-4.4, meanManuPair_pas_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot(-3.6, meanManuPair_pas_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot(-6.4, mean_pai_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot(-5.6, mean_pai_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot(-4.4, meanManuPair_pai_s_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot(-3.6, meanManuPair_pai_s_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')

    #reciprocal pairing, not threshold (a manual fascicle does not appear twice, nor an auto fascicle)
    axSD.plot(rand_jitter([-2-0.4]*len(d_s_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_s_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axSD.plot(rand_jitter([-2+0.4]*len(d_s_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_s_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axSD.plot([-2-0.4], meanUniqueMatch_d_s_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSD.plot([-2+0.4], meanUniqueMatch_d_s_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot(rand_jitter([-2-0.4]*len(fl_s_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_s_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axSfl.plot(rand_jitter([-2+0.4]*len(fl_s_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_s_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axSfl.plot([-2-0.4], meanUniqueMatch_fl_s_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSfl.plot([-2+0.4], meanUniqueMatch_fl_s_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot(rand_jitter([-2-0.4]*len(pas_s_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(pas_s_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axSpas.plot(rand_jitter([-2+0.4]*len(pas_s_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(pas_s_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axSpas.plot([-2-0.4], meanUniqueMatch_pas_s_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpas.plot([-2+0.4], meanUniqueMatch_pas_s_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot(rand_jitter([-2-0.4]*len(pai_s_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(pai_s_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axSpai.plot(rand_jitter([-2+0.4]*len(pai_s_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(pai_s_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axSpai.plot([-2-0.4], meanUniqueMatch_pai_s_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axSpai.plot([-2+0.4], meanUniqueMatch_pai_s_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')

    axSrel.plot(d_s_a_rel, fl_s_a_rel, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')
    axSrel.plot(d_s_m_rel, fl_s_m_rel, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')

    plt.show()
    figSD.savefig('C:/Users/Lisa Paillard/Desktop/testSDthresholds.jpg')
    figSrel.savefig('C:/Users/Lisa Paillard/Desktop/fl_DependingOn_d_s.jpg')
    figSfl.savefig('C:/Users/Lisa Paillard/Desktop/fl_dependingOn_thresholds_s.jpg')
    figSpas.savefig('C:/Users/Lisa Paillard/Desktop/pas_dependingOn_thresholds_s.jpg')
    figSpai.savefig('C:/Users/Lisa Paillard/Desktop/pai_dependingOn_thresholds_s.jpg')


def plotPairingPano(path_to_dict, name_dict, colors, participants):
    """Function that aims at comparing the way fascicles should be matched
    in panoramic images (in particular, maximum thresholds between match auto and manual fascicles)
    path_to_dict: path to folder containing the dictionnary
    name_dict: name of the dictionnary within the folder
    colors: list of color tuples; length should be at least equal to
    the number of analysed participants.
    participant: list of analysed participants' names/ID.     
    """    
    import numpy as np
    import fnmatch
    from dictmanager import load_obj
    
    #d = list of distances from reference point
    #fl = list of FL
    #PAs = list of PA sup
    #PAi = list of PA inf
    #mt = list of MT
    #_p : panoramic images
    #_m : manual
    #_a : automated

    '************************************************************************'
    '*****************************INITIALIZATION*****************************'
    
    spacial_thresh = [20,15,10,5,3,2,1,0] #millimeters    
    
    d_p_m = [[] for par in range(len(participants))]
    fl_p_m = [[] for par in range(len(participants))]
    PAs_p_m = [[] for par in range(len(participants))]
    PAi_p_m = [[] for par in range(len(participants))]
    mt_p_m = [[] for par in range(len(participants))]
    d_p_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    fl_p_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAs_p_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAi_p_m_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    
    d_p_a = [[] for par in range(len(participants))]
    fl_p_a = [[] for par in range(len(participants))]
    PAs_p_a = [[] for par in range(len(participants))]
    PAi_p_a = [[] for par in range(len(participants))]
    mt_p_a = [[] for par in range(len(participants))]
    d_p_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    fl_p_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAs_p_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]
    PAi_p_a_filtered = [[[] for par in range(len(participants))] for thr in range(len(spacial_thresh))]

    diff_calfct_p = []



    #plots to determine how to pair auto and manual fascicles
    manuPairing_d_p_m=[[] for par in range(len(participants))]
    manuPairing_d_p_a=[[] for par in range(len(participants))]
    manuPairing_fl_p_m=[[] for par in range(len(participants))]
    manuPairing_fl_p_a=[[] for par in range(len(participants))]
    manuPairing_pas_p_m=[[] for par in range(len(participants))]
    manuPairing_pas_p_a=[[] for par in range(len(participants))]
    manuPairing_pai_p_m=[[] for par in range(len(participants))]
    manuPairing_pai_p_a=[[] for par in range(len(participants))]
    d_p_m_rel=[]
    d_p_a_rel=[]
    fl_p_m_rel=[]
    fl_p_a_rel=[]
    pas_p_m_rel=[]
    pas_p_a_rel=[]
    pai_p_m_rel=[]
    pai_p_a_rel=[]
    
    
    
    figPD, axPD = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figPD.suptitle('Distance to the distal insertion point (mm) depending on the pairing, in panoramic images', fontsize=40)
    axPD.set_ylabel('Distance to distal insertion (in mm)', fontsize= 35)
    axPD.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axPD.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axPD.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axPD.grid(True)
    axPD.tick_params(axis='x', labelsize=28, colors='k')
    axPD.tick_params(axis='y', labelsize=25, colors='k')

    figPfl, axPfl = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figPfl.suptitle('FL (mm) depending on the pairing, in panoramic images', fontsize=40)
    axPfl.set_ylabel('FL (in mm)', fontsize= 35)
    axPfl.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axPfl.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axPfl.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axPfl.grid(True)
    axPfl.tick_params(axis='x', labelsize=28, colors='k')
    axPfl.tick_params(axis='y', labelsize=25, colors='k')
    
    figPrel, axPrel = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figPrel.suptitle('Influence of the distance from reference point on fascicle length estimation (panoramic images)', fontsize=40)
    axPrel.set_xlabel('Distance from reference point (in mm)', fontsize= 35)
    axPrel.set_ylabel('Manual and auto fascicles length (in mm)', fontsize= 35)
    axPrel.grid(True)
    axPrel.tick_params(axis='x', labelsize=28, colors='k')
    axPrel.tick_params(axis='y', labelsize=25, colors='k')

    figPpas, axPpas = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figPpas.suptitle('PA superior (degree) depending on the pairing, in panoramic images', fontsize=40)
    axPpas.set_ylabel('PA superior (degree)', fontsize= 35)
    axPpas.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axPpas.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axPpas.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axPpas.grid(True)
    axPpas.tick_params(axis='x', labelsize=28, colors='k')
    axPpas.tick_params(axis='y', labelsize=25, colors='k')

    figPpai, axPpai = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figPpai.suptitle('PA inferior (degree) depending on the pairing, in panoramic images', fontsize=40)
    axPpai.set_ylabel('PA inferior (degree)', fontsize= 35)
    axPpai.set_xlabel('Manual and auto fascicles per matching threshold', fontsize= 35)
    axPpai.set_xticks([-6,-4,-2,0,2,4,6,8,10,12,14])
    axPpai.set_xticklabels(['not paired', 'Best auto equivalent', 'reciprocal pairing','20mm', '15mm', '10mm', '5mm', '3mm', '2mm','1mm','0mm'], rotation=30, ha='right', fontsize=25)    
    axPpai.grid(True)
    axPpai.tick_params(axis='x', labelsize=28, colors='k')
    axPpai.tick_params(axis='y', labelsize=25, colors='k')

    '************************************************************************'
    '*****************************DATA RETRIEVAL*****************************'

    dictio = load_obj(name_dict, path_to_dict)
    l2 = ['fasc*', 'fsc_*']
    
    for par in range(len(participants)):
        
        participant = participants[par]
        fam_folders = [str(d) for d in dictio[participant].keys()]

        p_manuFasc = []
        p_autoFasc = []
        
        for fam in fam_folders:
            
            ###################################################################
            # panoramic images
            dictioP = dictio[participant][fam]['BF']['panoramic']
            images = [str(im) for im in dictioP.keys()]
            for i in images:
                
                ###############################################################
                # PANORAMIC - manual
                dictioM = dictioP[i]['architecture manual']
                fascicles = [fa for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    idf = fam + '/' + i + '/' + f
                    if len(dictioF.keys())>1:
                        p_manuFasc.append((idf, dictioF['dist from insertion in mm']))
                        d_p_m[par].append(dictioF['dist from insertion in mm'])
                        fl_p_m[par].append(dictioF['FL']['length in mm'])
                        PAs_p_m[par].append(dictioF['PAsup']['value in degree'])
                        PAi_p_m[par].append(dictioF['PAinf']['value in degree'])            

                ###############################################################
                # PANORAMIC - automatic
                if ('architecture auto' in dictioP[i]):
                    dictioA = dictioP[i]['architecture auto']
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        for f in fascicles:
                            dictioF = dictioA[f]
                            idf = fam + '/' + i + '/' + f
                            #only keep fascicles that are entirely within the cropped image,
                            #to compare with manually identified fascicles
                            if len(dictioF.keys())>1 and dictioF['FL']['in/out of the image'] == 'in image':
                                p_autoFasc.append((idf, dictioF['dist from insertion in mm']))
                                d_p_a[par].append(dictioF['dist from insertion in mm'])
                                fl_p_a[par].append(dictioF['FL']['length in mm'])
                                PAs_p_a[par].append(dictioF['PAsup']['value in degree'])
                                PAi_p_a[par].append(dictioF['PAinf']['value in degree'])
                                

                                
                        if ('MT for labelled points' in dictioM['MT']):
                            for ind0 in range(len(dictioM['MT']['MT for labelled points'])):
                                elem = dictioM['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_p_m[par].append(elem) #MT in mm
                            
                            for ind0 in range(len(dictioA['MT']['MT for labelled points'])):
                                elem = dictioA['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_p_a[par].append(elem)

                    #calibration factors difference
                    diff_calfct_p.append(abs(dictioM['calfct_to_mm']-dictioA['calfct_to_mm before resize']['horizontal axis']))

        '************************************************************************'
        '***************SEARCH FOR MATCHING AUTO & MANUAL FASCICLES**************'
        listePair_manuF_p = []

        for n in range(len(p_manuFasc)):
            mf = p_manuFasc[n]
            subtr = [(tup,abs(tup[1]- mf[1])) for tup in p_autoFasc]
            subtr.sort(key=lambda x:x[1])
            closest = subtr[0]
            listePair_manuF_p.append((mf[0], closest[0][0], closest[1])) #tuple = ( ID manu fasc, ID auto fasc, distance entre les deux)
        listePair_manuF_p.sort(key=lambda x:x[1])
        uniqueMatching = []
        counterL = 0
        while counterL < len(listePair_manuF_p):
            currentAutoFasc = listePair_manuF_p[counterL][1]
            correspondingAutoFasc = [(listePair_manuF_p[counterL][0], listePair_manuF_p[counterL][2])]
            rank = counterL + 1
            while rank<len(listePair_manuF_p) and listePair_manuF_p[rank][1] == currentAutoFasc:
                correspondingAutoFasc.append((listePair_manuF_p[rank][0],listePair_manuF_p[rank][2]))
                rank = rank + 1
            correspondingAutoFasc.sort(key=lambda x:x[1])
            uniqueMatching.append((correspondingAutoFasc[0][0], currentAutoFasc, correspondingAutoFasc[0][1]))
            counterL = rank

        
        #manual fascicles pairing with an auto fascicle (manual fascicles can be used for several auto fasc)  
        for el in listePair_manuF_p:
            pathA0 = el[1].split('/')
            pathM0 = el[0].split('/')
            manuPairing_d_p_m[par].append(dictio[participant][pathM0[0]]['BF']['panoramic'][pathM0[1]]['architecture manual'][pathM0[2]]['dist from insertion in mm'])
            manuPairing_d_p_a[par].append(dictio[participant][pathA0[0]]['BF']['panoramic'][pathA0[1]]['architecture auto'][pathA0[2]]['dist from insertion in mm'])
            manuPairing_fl_p_m[par].append(dictio[participant][pathM0[0]]['BF']['panoramic'][pathM0[1]]['architecture manual'][pathM0[2]]['FL']['length in mm'])
            manuPairing_fl_p_a[par].append(dictio[participant][pathA0[0]]['BF']['panoramic'][pathA0[1]]['architecture auto'][pathA0[2]]['FL']['length in mm'])
            manuPairing_pas_p_m[par].append(dictio[participant][pathM0[0]]['BF']['panoramic'][pathM0[1]]['architecture manual'][pathM0[2]]['PAsup']['value in degree'])
            manuPairing_pas_p_a[par].append(dictio[participant][pathA0[0]]['BF']['panoramic'][pathA0[1]]['architecture auto'][pathA0[2]]['PAsup']['value in degree'])
            manuPairing_pai_p_m[par].append(dictio[participant][pathM0[0]]['BF']['panoramic'][pathM0[1]]['architecture manual'][pathM0[2]]['PAinf']['value in degree'])
            manuPairing_pai_p_a[par].append(dictio[participant][pathA0[0]]['BF']['panoramic'][pathA0[1]]['architecture auto'][pathA0[2]]['PAinf']['value in degree'])
        
        #reciprocal matching without any threshold on the distance between auto and manual fascicles
        for element in uniqueMatching:
            pathA1 = element[1].split('/')
            pathM1 = element[0].split('/')
            d_p_m_rel.append(dictio[participant][pathM1[0]]['BF']['panoramic'][pathM1[1]]['architecture manual'][pathM1[2]]['dist from insertion in mm'])
            d_p_a_rel.append(dictio[participant][pathA1[0]]['BF']['panoramic'][pathA1[1]]['architecture auto'][pathA1[2]]['dist from insertion in mm'])
            fl_p_m_rel.append(dictio[participant][pathM1[0]]['BF']['panoramic'][pathM1[1]]['architecture manual'][pathM1[2]]['FL']['length in mm'])
            fl_p_a_rel.append(dictio[participant][pathA1[0]]['BF']['panoramic'][pathA1[1]]['architecture auto'][pathA1[2]]['FL']['length in mm'])
            pas_p_m_rel.append(dictio[participant][pathM1[0]]['BF']['panoramic'][pathM1[1]]['architecture manual'][pathM1[2]]['PAsup']['value in degree'])
            pas_p_a_rel.append(dictio[participant][pathA1[0]]['BF']['panoramic'][pathA1[1]]['architecture auto'][pathA1[2]]['PAsup']['value in degree'])
            pai_p_m_rel.append(dictio[participant][pathM1[0]]['BF']['panoramic'][pathM1[1]]['architecture manual'][pathM1[2]]['PAinf']['value in degree'])
            pai_p_a_rel.append(dictio[participant][pathA1[0]]['BF']['panoramic'][pathA1[1]]['architecture auto'][pathA1[2]]['PAinf']['value in degree'])
        
        
        #reciprocal pairing with threshold th
        for t in range(len(spacial_thresh)):
            th = spacial_thresh[t]
            for element in uniqueMatching:
                pathA = element[1].split('/')
                pathM = element[0].split('/')
                
                if element[2] < th :
    
                    d_p_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['dist from insertion in mm'])
                    fl_p_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['FL']['length in mm'])
                    PAs_p_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['PAsup']['value in degree'])
                    PAi_p_m_filtered[t][par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['PAinf']['value in degree'])
                    
                    d_p_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['dist from insertion in mm'])
                    fl_p_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['FL']['length in mm'])
                    PAs_p_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['PAsup']['value in degree'])
                    PAi_p_a_filtered[t][par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['PAinf']['value in degree'])

            #reciprocal pairing with threshold th
            axPD.plot(rand_jitter([2*t-0.4]*len(d_p_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_p_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axPD.plot(rand_jitter([2*t+0.4]*len(d_p_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_p_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axPfl.plot(rand_jitter([2*t-0.4]*len(fl_p_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_p_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axPfl.plot(rand_jitter([2*t+0.4]*len(fl_p_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_p_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axPpas.plot(rand_jitter([2*t-0.4]*len(PAs_p_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAs_p_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axPpas.plot(rand_jitter([2*t+0.4]*len(PAs_p_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAs_p_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
            axPpai.plot(rand_jitter([2*t-0.4]*len(PAi_p_m_filtered[t][par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAi_p_m_filtered[t][par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
            axPpai.plot(rand_jitter([2*t+0.4]*len(PAi_p_a_filtered[t][par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAi_p_a_filtered[t][par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

        #all manu fascicles were paired with the closest auto fasc -> one auto fasc could be used several time
        axPD.plot(rand_jitter([-4-0.4]*len(manuPairing_d_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_d_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPD.plot(rand_jitter([-4+0.4]*len(manuPairing_d_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_d_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPfl.plot(rand_jitter([-4-0.4]*len(manuPairing_fl_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_fl_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPfl.plot(rand_jitter([-4+0.4]*len(manuPairing_fl_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_fl_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPpas.plot(rand_jitter([-4-0.4]*len(manuPairing_pas_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_pas_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPpas.plot(rand_jitter([-4+0.4]*len(manuPairing_pas_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_pas_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPpai.plot(rand_jitter([-4-0.4]*len(manuPairing_pai_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(manuPairing_pai_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPpai.plot(rand_jitter([-4+0.4]*len(manuPairing_pai_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(manuPairing_pai_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

        #all fascicles, not even paired
        axPD.plot(rand_jitter([-6-0.4]*len(d_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPD.plot(rand_jitter([-6+0.4]*len(d_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPfl.plot(rand_jitter([-6-0.4]*len(fl_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPfl.plot(rand_jitter([-6+0.4]*len(fl_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPpas.plot(rand_jitter([-6-0.4]*len(PAs_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAs_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPpas.plot(rand_jitter([-6+0.4]*len(PAs_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAs_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
        axPpai.plot(rand_jitter([-6-0.4]*len(PAi_p_m[par]), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(PAi_p_m[par]), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
        axPpai.plot(rand_jitter([-6+0.4]*len(PAi_p_a[par]), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(PAi_p_a[par]), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   

    for t2 in range(len(spacial_thresh)):
        
        mean_d_p_a_th =np.mean([item for sublist in d_p_a_filtered[t2] for item in sublist])
        mean_d_p_m_th = np.mean([item for sublist in d_p_m_filtered[t2] for item in sublist])
        mean_fl_p_a_th =np.mean([item for sublist in fl_p_a_filtered[t2] for item in sublist])
        mean_fl_p_m_th = np.mean([item for sublist in fl_p_m_filtered[t2] for item in sublist])
        mean_pas_p_m_th=np.mean([item for sublist in PAs_p_m_filtered[t2] for item in sublist])
        mean_pas_p_a_th=np.mean([item for sublist in PAs_p_a_filtered[t2] for item in sublist])
        mean_pai_p_m_th=np.mean([item for sublist in PAi_p_m_filtered[t2] for item in sublist])
        mean_pai_p_a_th=np.mean([item for sublist in PAi_p_a_filtered[t2] for item in sublist])
        axPD.plot(2*t2-0.4, mean_d_p_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPD.plot(2*t2+0.4, mean_d_p_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPfl.plot(2*t2-0.4, mean_fl_p_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPfl.plot(2*t2+0.4, mean_fl_p_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPpas.plot(2*t2-0.4, mean_pas_p_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPpas.plot(2*t2+0.4, mean_pas_p_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPpai.plot(2*t2-0.4, mean_pai_p_m_th, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')   
        axPpai.plot(2*t2+0.4, mean_pai_p_a_th, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')   
    
    mean_d_p_m = np.mean([item for sublist in d_p_m for item in sublist])
    mean_d_p_a = np.mean([item for sublist in d_p_a for item in sublist])
    mean_fl_p_m = np.mean([item for sublist in fl_p_m for item in sublist])
    mean_fl_p_a = np.mean([item for sublist in fl_p_a for item in sublist])
    mean_pas_p_m = np.mean([item for sublist in PAs_p_m for item in sublist])
    mean_pas_p_a =np.mean([item for sublist in PAs_p_a for item in sublist])
    mean_pai_p_m=np.mean([item for sublist in PAi_p_m for item in sublist])
    mean_pai_p_a=np.mean([item for sublist in PAi_p_a for item in sublist])
    meanManuPair_d_p_m = np.mean([item for sublist in manuPairing_d_p_m for item in sublist])
    meanManuPair_d_p_a = np.mean([item for sublist in manuPairing_d_p_a for item in sublist])
    meanManuPair_fl_p_m = np.mean([item for sublist in manuPairing_fl_p_m for item in sublist])
    meanManuPair_fl_p_a = np.mean([item for sublist in manuPairing_fl_p_a for item in sublist])
    meanManuPair_pas_p_m = np.mean([item for sublist in manuPairing_pas_p_m for item in sublist])
    meanManuPair_pas_p_a = np.mean([item for sublist in manuPairing_pas_p_a for item in sublist])
    meanManuPair_pai_p_m = np.mean([item for sublist in manuPairing_pai_p_m for item in sublist])
    meanManuPair_pai_p_a = np.mean([item for sublist in manuPairing_pai_p_a for item in sublist])
    meanUniqueMatch_d_p_m = np.mean(d_p_m_rel)
    meanUniqueMatch_d_p_a = np.mean(d_p_a_rel)
    meanUniqueMatch_fl_p_m = np.mean(fl_p_m_rel)
    meanUniqueMatch_fl_p_a = np.mean(fl_p_a_rel)
    meanUniqueMatch_pas_p_m = np.mean(pas_p_m_rel)
    meanUniqueMatch_pas_p_a = np.mean(pas_p_a_rel)
    meanUniqueMatch_pai_p_m = np.mean(pai_p_m_rel)
    meanUniqueMatch_pai_p_a = np.mean(pai_p_a_rel)

    axPD.plot(-6.4, mean_d_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPD.plot(-5.6, mean_d_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPD.plot(-4.4, meanManuPair_d_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPD.plot(-3.6, meanManuPair_d_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot(-6.4, mean_fl_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot(-5.6, mean_fl_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot(-4.4, meanManuPair_fl_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot(-3.6, meanManuPair_fl_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot(-6.4, mean_pas_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot(-5.6, mean_pas_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot(-4.4, meanManuPair_pas_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot(-3.6, meanManuPair_pas_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot(-6.4, mean_pai_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot(-5.6, mean_pai_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot(-4.4, meanManuPair_pai_p_m, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot(-3.6, meanManuPair_pai_p_a, color = 'k', marker='x', markeredgewidth = 2, markersize = 60, linestyle='None')

    #reciprocal pairing, not threshold (a manual fascicle does not appear twice, nor an auto fascicle)
    axPD.plot(rand_jitter([-2-0.4]*len(d_p_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(d_p_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axPD.plot(rand_jitter([-2+0.4]*len(d_p_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(d_p_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axPD.plot([-2-0.4], meanUniqueMatch_d_p_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPD.plot([-2+0.4], meanUniqueMatch_d_p_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot(rand_jitter([-2-0.4]*len(fl_p_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(fl_p_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axPfl.plot(rand_jitter([-2+0.4]*len(fl_p_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(fl_p_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axPfl.plot([-2-0.4], meanUniqueMatch_fl_p_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPfl.plot([-2+0.4], meanUniqueMatch_fl_p_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot(rand_jitter([-2-0.4]*len(pas_p_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(pas_p_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axPpas.plot(rand_jitter([-2+0.4]*len(pas_p_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(pas_p_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axPpas.plot([-2-0.4], meanUniqueMatch_pas_p_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpas.plot([-2+0.4], meanUniqueMatch_pas_p_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot(rand_jitter([-2-0.4]*len(pai_p_m_rel), sensib=0.07, lowerLimit=-1, upperLimit=0.2), rand_jitter(pai_p_m_rel), alpha = 0.2, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')   
    axPpai.plot(rand_jitter([-2+0.4]*len(pai_p_a_rel), sensib=0.07, lowerLimit=-0.2, upperLimit=1), rand_jitter(pai_p_a_rel), alpha = 0.2, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')   
    axPpai.plot([-2-0.4], meanUniqueMatch_pai_p_m, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')
    axPpai.plot([-2+0.4], meanUniqueMatch_pai_p_a, color = 'k', marker='x',  markeredgewidth = 2, markersize = 60, linestyle='None')

    axPrel.plot(d_p_a_rel, fl_p_a_rel, color = (1,0.26,0.11), marker='.', markersize = 40, linestyle='None')
    axPrel.plot(d_p_m_rel, fl_p_m_rel, color = (0.2,0.8,0.2), marker='.', markersize = 40, linestyle='None')

    plt.show()
    figPD.savefig('C:/Users/Lisa Paillard/Desktop/testPDthresholds.jpg')
    figPrel.savefig('C:/Users/Lisa Paillard/Desktop/fl_DependingOn_d_p.jpg')
    figPfl.savefig('C:/Users/Lisa Paillard/Desktop/fl_dependingOn_thresholdp_p.jpg')
    figPpas.savefig('C:/Users/Lisa Paillard/Desktop/pas_dependingOn_thresholdp_p.jpg')
    figPpai.savefig('C:/Users/Lisa Paillard/Desktop/pai_dependingOn_thresholdp_p.jpg')


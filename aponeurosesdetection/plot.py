"""plot"""
#colors
col = [(1,0.26,0.11),(0.45,0.45,1),(0.8,0.2,1),(0,0.4,1),\
          (0,1,0.6),(0.2,0.8,0.2),(1,0.8,0),(0.8,0,0),\
          (1,0,0.4),(0.7,0.24,0),(0,0,0.8)]

def plotFeatures(path_to_dict, name_dict, colors, participants):
    import numpy as np
    import fnmatch
    from dictmanager import load_obj
    
    
    #d = list of distances from reference point
    #fl = list of FL
    #PAs = list of PA sup
    #PAi = list of PA inf
    #mt = list of MT
    #_s : simple images
    #_p : panoramic images
    #_m : manual
    #_a : automated
    
    d_s_m = [[] for par in range(len(participants))]
    fl_s_m = [[] for par in range(len(participants))]
    PAs_s_m = [[] for par in range(len(participants))]
    PAi_s_m = [[] for par in range(len(participants))]
    mt_s_m = [[] for par in range(len(participants))]
    
    d_s_a = [[] for par in range(len(participants))]
    fl_s_a = [[] for par in range(len(participants))]
    PAs_s_a = [[] for par in range(len(participants))]
    PAi_s_a = [[] for par in range(len(participants))]
    mt_s_a = [[] for par in range(len(participants))]
    
    d_p_m = [[] for par in range(len(participants))]
    fl_p_m = [[] for par in range(len(participants))]
    PAs_p_m = [[] for par in range(len(participants))]
    PAi_p_m = [[] for par in range(len(participants))]
    mt_p_m = [[] for par in range(len(participants))]
    
    d_p_a = [[] for par in range(len(participants))]
    fl_p_a = [[] for par in range(len(participants))]
    PAs_p_a = [[] for par in range(len(participants))]
    PAi_p_a = [[] for par in range(len(participants))]
    mt_p_a = [[] for par in range(len(participants))]

    diff_calfct_p = []
    diff_calfct_s = []
    mt_s_diff = []
    mt_p_diff = []
    mt_s_relative_diff = []
    mt_p_relative_diff = []

    dictio = load_obj(name_dict, path_to_dict)
    l2 = ['fasc*', 'fsc_*']
    
    for par in range(len(participants)):
        
        participant = participants[par]
        fam_folders = [str(d) for d in dictio[participant].keys()]

        for fam in fam_folders:
            # simple images
            dictioS = dictio[participant][fam]['BF']['simple']
            images = [str(im) for im in dictioS.keys()]
            for i in images:
                
                # manual
                dictioM = dictioS[i]['architecture manual']
                fascicles = [str(fa) for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    if len(dictioF.keys())>1:
                        d_s_m[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                        fl_s_m[par].append(dictioF['FL']['length in mm'])
                        PAs_s_m[par].append(dictioF['PAsup']['value in degree'])
                        PAi_s_m[par].append(dictioF['PAinf']['value in degree'])

                # automatic
                if ('architecture auto' in dictioS[i]):
                    dictioA = dictioS[i]['architecture auto']
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        for f in fascicles:
                            dictioF = dictioA[f]
                            if len(dictioF.keys())>1:
                                d_s_a[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                                fl_s_a[par].append(dictioF['FL']['length in mm'])
                                PAs_s_a[par].append(dictioF['PAsup']['value in degree'])
                                PAi_s_a[par].append(dictioF['PAinf']['value in degree'])
                        if ('MT for labelled points' in dictioM['MT']):
                            mt_s_m[par].append(dictioM['MT']['MT for labelled points']) #(absc,MT) in mm
                            mt_s_a[par].append(dictioA['MT']['MT for labelled points'])
                    
                    #calibration factors difference
                    diff_calfct_s.append(abs(dictioM['calfct_to_mm']-\
                                             dictioA['calfct_to_mm']['horizontal axis']))

            # panoramic images
            dictioP = dictio[participant][fam]['BF']['panoramic']
            images = [str(im) for im in dictioP.keys()]
            for i in images:
        
                # manual
                dictioM = dictioP[i]['architecture manual']
                fascicles = [fa for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    if len(dictioF.keys())>1:
                        d_p_m[par].append(dictioF['dist from insertion in mm'])
                        fl_p_m[par].append(dictioF['FL']['length in mm'])
                        PAs_p_m[par].append(dictioF['PAsup']['value in degree'])
                        PAi_p_m[par].append(dictioF['PAinf']['value in degree'])            
    
                # automatic
                if ('architecture auto' in dictioP[i]):
                    dictioA = dictioP[i]['architecture auto']
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        for f in fascicles:
                            dictioF = dictioA[f]
                            if len(dictioF.keys())>1 and dictioF['FL']['in/out of the image'] == 'in image':
                                d_p_a[par].append(dictioF['dist from insertion in mm'])
                                fl_p_a[par].append(dictioF['FL']['length in mm'])
                                PAs_p_a[par].append(dictioF['PAsup']['value in degree'])
                                PAi_p_a[par].append(dictioF['PAinf']['value in degree'])      
                        if ('MT for labelled points' in dictioM['MT']):
                            mt_p_m[par].append(dictioM['MT']['MT for labelled points'])
                            mt_p_a[par].append(dictioA['MT']['MT for labelled points'])
    
                    #calibration factors difference
                    diff_calfct_p.append(abs(dictioM['calfct_to_mm']-dictioA['calfct_to_mm before resize']['horizontal axis']))
                
    import matplotlib.pyplot as plt
    figS1, (axS1a, axS1b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figS2, (axS2a, axS2b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figS3, (axS3a, axS3b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figS4, axS4 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))

    figP1, (axP1a, axP1b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figP2, (axP2a, axP2b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figP3, (axP3a, axP3b) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(25,25))
    figP4, axP4 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))

    figCs, (ax1Cs, ax2Cs) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15))
    figCp, (ax1Cp, ax2Cp) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15))

    #muscle features
    figS1.suptitle('PA sup in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS1a.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axS1a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axS1a.grid(True)
    axS1b.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axS1b.grid(True)
    
    figS2.suptitle('PA inf in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS2a.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axS2a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axS2a.grid(True)
    axS2b.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axS2b.grid(True)
    
    figS3.suptitle('FL in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS3a.set_ylabel('FL (mm)', fontsize= 8)
    axS3a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axS3a.grid(True)
    axS3b.set_ylabel('FL (mm)', fontsize= 8)
    axS3b.grid(True)

    figS4.suptitle('MT in simple images')
    axS4.set_ylabel('MT (mm)', fontsize= 8)
    axS4.grid(True)

    figP1.suptitle('PA sup in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP1a.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axP1a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axP1a.grid(True)
    axP1b.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axP1b.grid(True)

    figP2.suptitle('PA inf in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP2a.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axP2a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axP2a.grid(True)
    axP2b.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axP2b.grid(True)

    figP3.suptitle('FL in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP3a.set_ylabel('FL (mm)', fontsize= 8)
    axP3a.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    axP3a.grid(True)
    axP3b.set_ylabel('FL (mm)', fontsize= 8)
    axP3b.grid(True)

    figP4.suptitle('MT in panoramic images')
    axP4.set_ylabel('MT (mm)', fontsize= 8)
    axP4.grid(True)

    labels = ['TOT-manu']
    for par in range(len(participants)):

        #statistics per participant
        median_p_a_PAsup = np.median(PAs_p_a[par])
        median_p_a_PAinf = np.median(PAi_p_a[par])
        median_p_a_FL = np.median(fl_p_a[par])
        median_p_m_PAsup = np.median(PAs_p_m[par])
        median_p_m_PAinf = np.median(PAi_p_m[par])
        median_p_m_FL = np.median(fl_p_m[par])
        median_s_a_PAsup =np.median(PAs_s_a[par])
        median_s_a_PAinf =np.median(PAi_s_a[par])
        median_s_a_FL =np.median(fl_s_a[par])
        median_s_m_PAsup =np.median(PAs_s_m[par])
        median_s_m_PAinf =np.median(PAi_s_m[par])
        median_s_m_FL =np.median(fl_s_m[par])
        
        print('part = ', par, 'stat=', median_p_a_PAsup, median_p_a_PAinf, median_p_a_FL, median_p_m_PAsup,\
              median_p_m_PAinf, median_p_m_FL, median_s_a_PAsup, median_s_a_PAinf,\
              median_s_a_FL, median_s_m_PAsup, median_s_m_PAinf, median_s_m_FL)
        
        #plots
        ID = participants[par][:2]
        labels.append(ID+'-auto')
        labels.append(ID+'-manu')
        
        axS3a.plot(d_s_a[par], fl_s_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axS3a.plot(d_s_m[par], fl_s_m[par], color = colors[par], marker='+', markersize = 7, linestyle='None')
        axS3b.plot(['TOT-auto']*len(fl_s_a[par]), fl_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS3b.plot(['TOT-manu']*len(fl_s_m[par]), fl_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS3b.plot([ID+'-auto']*len(fl_s_a[par]), fl_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS3b.plot([ID+'-manu']*len(fl_s_m[par]), fl_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS3b.plot([ID+'-auto'], median_s_a_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
        axS3b.plot([ID+'-manu'], median_s_m_FL, color = 'k', marker='x', markersize = 12, linestyle='None')

        axS1a.plot(d_s_a[par], PAs_s_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axS1a.plot(d_s_m[par], PAs_s_m[par], color = colors[par], marker='+', markersize = 7, linestyle='None')
        axS1b.plot(['TOT-auto']*len(PAs_s_a[par]), PAs_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS1b.plot(['TOT-manu']*len(PAs_s_m[par]), PAs_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS1b.plot([ID+'-auto']*len(PAs_s_a[par]), PAs_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS1b.plot([ID+'-manu']*len(PAs_s_m[par]), PAs_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS1b.plot([ID+'-auto'], median_s_a_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
        axS1b.plot([ID+'-manu'], median_s_m_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')

        axS2a.plot(d_s_a[par], PAi_s_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axS2a.plot(d_s_m[par], PAi_s_m[par], color = colors[par], marker='+', markersize = 7, linestyle='None')
        axS2b.plot(['TOT-auto']*len(PAi_s_a[par]), PAi_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS2b.plot(['TOT-manu']*len(PAi_s_m[par]), PAi_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS2b.plot([ID+'-auto']*len(PAi_s_a[par]), PAi_s_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS2b.plot([ID+'-manu']*len(PAi_s_m[par]), PAi_s_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axS2b.plot([ID+'-auto'], median_s_a_PAinf, 'k', marker='x', markersize = 12, linestyle='None')
        axS2b.plot([ID+'-manu'], median_s_m_PAinf, 'k', marker='x', markersize = 12, linestyle='None')

        for ind0 in range(len(mt_s_m[par])):
            for ind1 in range(len(mt_s_m[par][ind0])):
                elem2 = mt_s_m[par][ind0][ind1]
                axS4.plot(['manual'], [elem2], color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 0.5)
                if len(mt_s_a[par]) >0:
                    elem1 = mt_s_a[par][ind0][ind1]
                    if elem1 != 'error' and elem2 != 'error':
                        mt_s_diff.append(elem1-elem2)
                        mt_s_relative_diff.append(abs(elem1-elem2)/elem2)
                        axS4.plot(['automatic', 'manual'], [elem1, elem2], color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 0.5)
                    if 'TOT-auto' not in labels:
                        labels.append('TOT-auto')

        axP3a.plot(d_p_a[par], fl_p_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axP3a.plot(d_p_m[par], fl_p_m[par],  color = colors[par], marker='+', markersize = 7, linestyle='None')
        axP3b.plot(['TOT-auto']*len(fl_p_a[par]), fl_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP3b.plot(['TOT-manu']*len(fl_p_m[par]), fl_p_m[par],  color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP3b.plot([ID+'-auto']*len(fl_p_a[par]), fl_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP3b.plot([ID+'-manu']*len(fl_p_m[par]), fl_p_m[par],  color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP3b.plot([ID+'-auto'], median_p_a_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
        axP3b.plot([ID+'-manu'], median_p_m_FL, color = 'k', marker='x', markersize = 12, linestyle='None')

        axP1a.plot(d_p_a[par], PAs_p_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axP1a.plot(d_p_m[par], PAs_p_m[par], color = colors[par], marker='+', markersize = 7, linestyle='None')
        axP1b.plot(['TOT-auto']*len(PAs_p_a[par]), PAs_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP1b.plot(['TOT-manu']*len(PAs_p_m[par]), PAs_p_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP1b.plot([ID+'-auto']*len(PAs_p_a[par]), PAs_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP1b.plot([ID+'-manu']*len(PAs_p_m[par]), PAs_p_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP1b.plot([ID+'-auto'], median_p_a_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
        axP1b.plot([ID+'-manu'], median_p_m_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')

        axP2a.plot(d_p_a[par], PAi_p_a[par], color = colors[par], marker='.', markersize = 5, linestyle='None')
        axP2a.plot(d_p_m[par], PAi_p_m[par], color = colors[par], marker='+', markersize = 7, linestyle='None')
        axP2b.plot(['TOT-auto']*len(PAi_p_a[par]), PAi_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP2b.plot(['TOT-manu']*len(PAi_p_m[par]), PAi_p_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP2b.plot([ID+'-auto']*len(PAi_p_a[par]), PAi_p_a[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP2b.plot([ID+'-manu']*len(PAi_p_m[par]), PAi_p_m[par], color = colors[par], marker='.', markersize = 7, linestyle='None')
        axP2b.plot([ID+'-auto'], median_p_a_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')
        axP2b.plot([ID+'-manu'], median_p_m_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')

        for ind0 in range(len(mt_p_m[par])):
            for ind1 in range(len(mt_p_m[par][ind0])):
                elem2 = mt_p_m[par][ind0][ind1]
                axP4.plot(['manual'], [elem2], color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 0.5)
                if len(mt_p_a[par])>0:
                    elem1 = mt_p_a[par][ind0][ind1]
                    if elem1 != 'error' and elem2 != 'error':
                        mt_p_diff.append(elem1-elem2)
                        mt_p_relative_diff.append(abs(elem1-elem2)/elem2)
                        axP4.plot(['automatic', 'manual'], [elem1, elem2], color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 0.5)
                    if 'TOT-auto' not in labels:
                        labels.append('TOT-auto')
                        
    #statistics on the population
    median_p_a_PAsup = np.median([item for sublist in PAs_p_a for item in sublist])
    median_p_a_PAinf = np.median([item for sublist in PAi_p_a for item in sublist])
    median_p_a_FL = np.median([item for sublist in fl_p_a for item in sublist])
    median_p_m_PAsup = np.median([item for sublist in PAs_p_m for item in sublist])
    median_p_m_PAinf = np.median([item for sublist in PAi_p_m for item in sublist])
    median_p_m_FL = np.median([item for sublist in fl_p_m for item in sublist])
    median_s_a_PAsup =np.median([item for sublist in PAs_s_a for item in sublist])
    median_s_a_PAinf =np.median([item for sublist in PAi_s_a for item in sublist])
    median_s_a_FL =np.median([item for sublist in fl_s_a for item in sublist])
    median_s_m_PAsup =np.median([item for sublist in PAs_s_m for item in sublist])
    median_s_m_PAinf =np.median([item for sublist in PAi_s_m for item in sublist])
    median_s_m_FL =np.median([item for sublist in fl_s_m for item in sublist])

    axS3b.plot(['TOT-auto'], median_s_a_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
    axS3b.plot(['TOT-manu'], median_s_m_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
    axS1b.plot(['TOT-auto'], median_s_a_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
    axS1b.plot(['TOT-manu'], median_s_m_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
    axS2b.plot(['TOT-auto'], median_s_a_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')
    axS2b.plot(['TOT-manu'], median_s_m_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP3b.plot(['TOT-auto'], median_p_a_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP3b.plot(['TOT-manu'], median_p_m_FL, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP1b.plot(['TOT-auto'], median_p_a_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP1b.plot(['TOT-manu'], median_p_m_PAsup, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP2b.plot(['TOT-auto'], median_p_a_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')
    axP2b.plot(['TOT-manu'], median_p_m_PAinf, color = 'k', marker='x', markersize = 12, linestyle='None')
    
    print('all participants stats=', median_p_a_PAsup, median_p_a_PAinf, median_p_a_FL, median_p_m_PAsup,\
          median_p_m_PAinf, median_p_m_FL, median_s_a_PAsup, median_s_a_PAinf,\
          median_s_a_FL, median_s_m_PAsup, median_s_m_PAinf, median_s_m_FL)

    if len(mt_s_diff)>0 and len(mt_s_relative_diff)>0:
        mean_diff_mt_s = np.mean(mt_s_diff)
        mean_relative_diff_mt_s = np.mean(mt_s_relative_diff)
        max_relative_diff_mt_s = np.amax(mt_s_relative_diff)
        max_abs_diff_mt_s = np.amax(np.abs(mt_s_diff))
        median_diff_mt_s = np.median(mt_s_diff)
        print('Simple - Median != in MT (mm):', median_diff_mt_s)
        print('Simple - Mean != in MT (mm):', mean_diff_mt_s)
        print('Simple - Max absolute != in MT (mm):', max_abs_diff_mt_s)
        print('Simple - Mean %!= in MT:', mean_relative_diff_mt_s*100, ' %')
        print('Simple - Max %!= in MT:', max_relative_diff_mt_s*100, ' %')
        
    if len(mt_p_diff)>0 and len(mt_p_relative_diff)>0:
        mean_diff_mt_p = np.mean(mt_p_diff)
        mean_relative_diff_mt_p = np.mean(mt_p_relative_diff)
        max_abs_diff_mt_p = np.amax(np.abs(mt_p_diff))
        max_relative_diff_mt_p = np.amax(mt_p_relative_diff)
        median_diff_mt_p = np.median(mt_p_diff)
        print('Pano - Median != in MT (mm):', median_diff_mt_p)
        print('Pano - Mean != in MT (mm):', mean_diff_mt_p)
        print('Pano - Max absolute != in MT (mm):', max_abs_diff_mt_p)
        print('Pano - Mean %!= in MT:', mean_relative_diff_mt_p*100, ' %')
        print('Pano - Max %!= in MT:', max_relative_diff_mt_p*100, ' %')

    axS1b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    axS2b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    axS3b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    axP1b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    axP2b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)   
    axP3b.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)

    #calibration factors
    if len(diff_calfct_s)>0:
        diff_min_s = min(diff_calfct_s)
        diff_max_s = max(diff_calfct_s)
        diff_mean_s = np.mean(diff_calfct_s)
        figCs.suptitle('Difference in manual/automatic calibration factor in simple images')
        ax1Cs.plot([i for i in range(len(diff_calfct_s))], diff_calfct_s, color = colors[1], marker = 'o', markersize = 5, linestyle = 'None')
        ax1Cs.set_ylabel('Difference of calibration factors (mm/pixel)', fontsize= 8)
        ax1Cs.set_xlabel('image', fontsize= 8)
        ax1Cs.grid(True)
        ax2Cs.text(0.2,0.2, 'Mean ='+str(diff_mean_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cs.text(0.2,0.7, 'Min ='+str(diff_min_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cs.text(0.2,0.5, 'Max ='+str(diff_max_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cs.set_axis_off()

    if len(diff_calfct_p) >0:
        diff_min_p = min(diff_calfct_p)
        diff_max_p = max(diff_calfct_p)
        diff_mean_p = np.mean(diff_calfct_p)
        figCp.suptitle('Difference in manual/automatic calibration factor in panoramic images')
        ax1Cp.plot([i for i in range(len(diff_calfct_p))], diff_calfct_p, color = colors[1], marker = 'o', markersize = 5, linestyle = 'None')
        ax1Cp.set_ylabel('Difference of calibration factors (mm/pixel)', fontsize= 8)
        ax1Cp.set_xlabel('image', fontsize= 8)
        ax1Cp.grid(True)
        ax2Cp.text(0.2,0.2, 'Mean ='+str(diff_mean_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cp.text(0.2,0.7, 'Min ='+str(diff_min_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cp.text(0.2,0.5, 'Max ='+str(diff_max_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
        ax2Cp.set_axis_off()

        
    plt.show()
"""plot"""


def plotFeatures(path_to_dict, name_dict):

    import fnmatch
    from dictmanager import load_obj
    
    participants = ['01_Kevin', '02_rafaelopes']
    
    #d = list of distances from reference point
    #fl = list of FL
    #z1 = list of PA sup
    #z2 = list of PA inf
    #mt = list of arrays (absc, MT). One array per image
    #_s : simple images
    #_p : panoramic images
    #_m : manual
    #_a : automated
    
    d_s_m = [[] for par in range(len(participants))]
    fl_s_m = [[] for par in range(len(participants))]
    z1_s_m = [[] for par in range(len(participants))]
    z2_s_m = [[] for par in range(len(participants))]
    mt_s_m = [[] for par in range(len(participants))]
    
    d_s_a = [[] for par in range(len(participants))]
    fl_s_a = [[] for par in range(len(participants))]
    z1_s_a = [[] for par in range(len(participants))]
    z2_s_a = [[] for par in range(len(participants))]
    mt_s_a = [[] for par in range(len(participants))]
    
    d_p_m = [[] for par in range(len(participants))]
    fl_p_m = [[] for par in range(len(participants))]
    z1_p_m = [[] for par in range(len(participants))]
    z2_p_m = [[] for par in range(len(participants))]
    mt_p_m = [[] for par in range(len(participants))]
    
    d_p_a = [[] for par in range(len(participants))]
    fl_p_a = [[] for par in range(len(participants))]
    z1_p_a = [[] for par in range(len(participants))]
    z2_p_a = [[] for par in range(len(participants))]
    mt_p_a = [[] for par in range(len(participants))]

    diff_calfct_p = []
    diff_calfct_s = []

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
                c_start = dictioM['MT']['columns interval'][0] * dictioM['calfct_to_mm']
                c_end = dictioM['MT']['columns interval'][1] * dictioM['calfct_to_mm']
                mt_s_m[par].append(dictioM['MT']['coords']) #(absc,MT) in mm
                fascicles = [str(fa) for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    d_s_m[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                    fl_s_m[par].append(dictioF['FL']['length in mm'])
                    z1_s_m[par].append(dictioF['PAsup']['value in degree'])
                    z2_s_m[par].append(dictioF['PAinf']['value in degree'])

                    
                # automatic
                dictioA = dictioS[i]['architecture auto']
                if dictioA:
                    fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                    for f in fascicles:
                        dictioF = dictioA[f]
                        d_s_a[par].append(dictioF['dist from (0,0) of RGB image, in mm'])
                        fl_s_a[par].append(dictioF['FL']['length in mm'])
                        z1_s_a[par].append(dictioF['PAsup']['value in degree'])
                        z2_s_a[par].append(dictioF['PAinf']['value in degree'])

                 
                    #extract MT on the same interval than for manual processing
                    c_start_p = dictioA['crop']['columns'][0] +\
                                dictioA['MT']['columns interval'][0]                   
                    
                    dictioA['MT']['coords'][:,0] = dictioA['MT']['coords'][:,0]+\
                            c_start_p * dictioA['calfct_to_mm']['horizontal axis']
                            
                    for ind in range(dictioA['MT']['coords'].shape[0]):
                        if dictioA['MT']['coords'][ind,0] >= c_start and dictioA['MT']['coords'][ind,0]<= c_end:
                            mt_s_a.append(dictioA['MT']['coords'][ind,:])
                
                #calibration factors difference
                diff_calfct_s.append(abs(dictioM['calfct_to_mm']-\
                                         dictioA['calfct_to_mm']['horizontal axis']))

            # panoramic images
            dictioP = dictio[participant][fam]['BF']['panoramic']
            images = [str(im) for im in dictioP.keys()]
            for i in images:
        
                # manual
                dictioM = dictioP[i]['architecture manual']
                c_start = dictioM['MT']['columns interval'][0] * dictioM['calfct_to_mm']
                c_end = dictioM['MT']['columns interval'][1] * dictioM['calfct_to_mm']
                mt_p_m[par].append(dictioM['MT']['coords'])#(absc,MT) in mm             
                fascicles = [fa for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioM[f]
                    d_p_m[par].append(dictioF['dist from insertion in mm'])
                    fl_p_m[par].append(dictioF['FL']['length in mm'])
                    z1_p_m[par].append(dictioF['PAsup']['value in degree'])
                    z2_p_m[par].append(dictioF['PAinf']['value in degree'])            
    
                # automatic
                dictioA = dictioP[i]['architecture auto']
                if dictioA and ('MT' in dictioA):
                    fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                    for f in fascicles:
                        dictioF = dictioA[f]
                        d_p_a[par].append(dictioF['dist from insertion in mm'])
                        fl_p_a[par].append(dictioF['FL']['length in mm'])
                        z1_p_a[par].append(dictioF['PAsup']['value in degree'])
                        z2_p_a[par].append(dictioF['PAinf']['value in degree'])      
                   
                    #extract MT on the same interval than for manual processing
                    c_start_p = dictioA['crop']['columns'][0] +\
                                dictioA['MT']['columns interval'][0]                   
                    
                    dictioA['MT']['coords'][:,0] = dictioA['MT']['coords'][:,0]+\
                            c_start_p * dictioA['calfct_to_mm after resize']['horizontal axis']
                            
                    for ind in range(dictioA['MT']['coords'].shape[0]):
                        if dictioA['MT']['coords'][ind,0] >= c_start and dictioA['MT']['coords'][ind,0]<= c_end:
                            mt_p_a.append(dictioA['MT']['coords'][ind,:])

                #calibration factors difference
                diff_calfct_p.append(abs(dictioM['calfct_to_mm']-dictioA['calfct_to_mm before resize']['horizontal axis']))
                
                
    
    import matplotlib.pyplot as plt
    figS1, axS1 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figS2, axS2 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figS3, axS3 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figS4, axS4 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))

    figP1, axP1 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figP2, axP2 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figP3, axP3 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    figP4, axP4 = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))

    figCs, (ax1Cs, ax2Cs) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15))
    figCp, (ax1Cp, ax2Cp) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15))

    #calibration factors
    import numpy as np
    diff_min_s = min(diff_calfct_s)
    diff_max_s = max(diff_calfct_s)
    diff_mean_s = np.mean(diff_calfct_s)
    diff_min_p = min(diff_calfct_p)
    diff_max_p = max(diff_calfct_p)
    diff_mean_p = np.mean(diff_calfct_p)
    
    figCs.suptitle('Difference in manual/automatic calibration factor estimation in simple images')
    ax1Cs.plot([i for i in range(len(diff_calfct_p))], diff_calfct_p, 'bo', markersize = 5)
    ax1Cs.set_ylabel('Difference of calibration factors (mm/pixel)', fontsize= 8)
    ax1Cs.set_xlabel('image', fontsize= 8)
    ax2Cs.text(0.2,0.2, 'Mean ='+str(diff_mean_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cs.text(0.2,0.7, 'Min ='+str(diff_min_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cs.text(0.2,0.5, 'Max ='+str(diff_max_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cs.set_axis_off()

    figCp.suptitle('Difference in manual/automatic calibration factor estimation in panoramic images')
    ax1Cp.plot([i for i in range(len(diff_calfct_p))], diff_calfct_p, 'bo', markersize = 5)
    ax1Cp.set_ylabel('Difference of calibration factors (mm/pixel)', fontsize= 8)
    ax1Cp.set_xlabel('image', fontsize= 8)
    ax2Cp.text(0.2,0.2, 'Mean ='+str(diff_mean_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cp.text(0.2,0.7, 'Min ='+str(diff_min_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cp.text(0.2,0.5, 'Max ='+str(diff_max_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
    ax2Cp.set_axis_off()

    
    #muscle features
    figS1.suptitle('PA sup in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS1.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axS1.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figS2.suptitle('PA inf in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS2.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axS2.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figS3.suptitle('FL in simple images. Points: auto data, Plus: manual data. One color/participant')
    axS3.set_ylabel('FL (mm)', fontsize= 8)
    axS3.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figS4.suptitle('MT in simple images')
    axS4.set_ylabel('MT (mm)', fontsize= 8)
    axS4.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figP1.suptitle('PA sup in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP1.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    axP1.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figP2.suptitle('PA inf in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP2.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    axP2.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    figP3.suptitle('FL in panoramic images. Points: auto data, Plus: manual data. One color/participant')
    axP3.set_ylabel('FL (mm)', fontsize= 8)
    axP3.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    
    figP4.suptitle('MT in panoramic images')
    axP4.set_ylabel('MT (mm)', fontsize= 8)
    axP4.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)

    #colors
    colors = [(0,0,0), (1,1,0), (1,0,1), (0,1,1), (0,1,0), (1,0,0), (0,0,1),\
                 (0.5,0.8,0),(0.5,0.8,0.5), (0.1,0.8,0),(0.1,0.2,0.1), (1,0.2,0),\
                (0.3,0.3,1), (1,0.15,0.15), (0,0.2,0.9), (0,0.2,0.15)]

    for par in range(len(participants)):
          
        axS3.plot(d_s_a[par], fl_s_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axS3.plot(d_s_m[par], fl_s_m[par], color = colors[par], marker='+', markersize = 5, linestyle='None')

        axS1.plot(d_s_a[par], z1_s_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axS1.plot(d_s_m[par], z1_s_m[par], color = colors[par], marker='+', markersize = 5, linestyle='None')

        axS2.plot(d_s_a[par], z2_s_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axS2.plot(d_s_m[par], z2_s_m[par], color = colors[par], marker='+', markersize = 5, linestyle='None')

        axP3.plot(d_p_a[par], fl_p_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axP3.plot(d_p_m[par], fl_p_m[par],  color = colors[par], marker='+', markersize = 5, linestyle='None')
    
        axP1.plot(d_p_a[par], z1_p_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axP1.plot(d_p_m[par], z1_p_m[par], color = colors[par], marker='+', markersize = 5, linestyle='None')
    
        axP2.plot(d_p_a[par], z2_p_a[par], color = colors[par], marker='.', markersize = 3, linestyle='None')
        axP2.plot(d_p_m[par], z2_p_m[par], color = colors[par], marker='+', markersize = 5, linestyle='None')

    plt.show()   
"""plot"""


def plotFeatures(participant, path_to_dict, name_dict):

    import fnmatch
    from dictmanager import load_obj
    
    x_s_m = []
    y_s_m = []
    z1_s_m = []
    z2_s_m = []
    x_s_a = []
    y_s_a = []
    z1_s_a = []
    z2_s_a = []
    x_p_m = []
    y_p_m = []
    z1_p_m = []
    z2_p_m = []
    x_p_a = []
    y_p_a = []
    z1_p_a = []
    z2_p_a = []

    dictio = load_obj(name_dict, path_to_dict)
    l2 = ['fasc*', 'fsc_*']
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
                x_s_m.append(dictioF['dist from (0,0) of RGB image, in mm'])
                y_s_m.append(dictioF['FL']['length in mm'])
                z1_s_m.append(dictioF['PAsup']['value in degree'])
                z2_s_m.append(dictioF['PAinf']['value in degree'])
                
            # automatic
            dictioA = dictioS[i]['architecture auto']
            if dictioA:
                fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioA[f]
                    x_s_a.append(dictioF['dist from (0,0) of RGB image, in mm'])
                    y_s_a.append(dictioF['FL']['length in mm'])
                    z1_s_a.append(dictioF['PAsup']['value in degree'])
                    z2_s_a.append(dictioF['PAinf']['value in degree'])        
        
        # panoramic images
        dictioP = dictio[participant][fam]['BF']['panoramic']
        images = [str(im) for im in dictioP.keys()]
        for i in images:
    
            # manual
            dictioM = dictioP[i]['architecture manual']
            fascicles = [fa for fa in dictioM if any(fnmatch.fnmatch(fa, p) for p in l2)]
            for f in fascicles:
                dictioF = dictioM[f]
                x_p_m.append(dictioF['dist from insertion in mm'])
                y_p_m.append(dictioF['FL']['length in mm'])
                z1_p_m.append(dictioF['PAsup']['value in degree'])
                z2_p_m.append(dictioF['PAinf']['value in degree'])            
            print('PM', len(x_p_m), len(y_p_m), len(z1_p_m), len(z2_p_m))

            # automatic
            dictioA = dictioP[i]['architecture auto']
            if dictioA:
                fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                for f in fascicles:
                    dictioF = dictioA[f]
                    x_p_a.append(dictioF['dist from insertion in mm'])
                    y_p_a.append(dictioF['FL']['length in mm'])
                    z1_p_a.append(dictioF['PAsup']['value in degree'])
                    z2_p_a.append(dictioF['PAinf']['value in degree'])        
                print('PA', len(x_p_a), len(y_p_a), len(z1_p_a), len(z2_p_a))

    import matplotlib.pyplot as plt
    figS, ((ax1S, ax2S), (ax3S, ax4S)) = plt.subplots(2,2,sharex = False, sharey = False, figsize =(15,15))
    figP, ((ax1P, ax2P), (ax3P, ax4P)) = plt.subplots(2,2,sharex = False, sharey = False, figsize =(15,15))
    
    #simple
    figS.suptitle('Simple images outputs for participant: ' + participant)
    
    ax1S.plot(x_s_a, y_s_a, 'b+', markersize = 3, label = 'auto')
    ax1S.plot(x_s_m, y_s_m, 'r+', markersize = 3, label = 'manual')
    ax1S.set_ylabel('FL (mm)', fontsize= 8)
    ax1S.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    ax1S.legend(loc = 'upper left', prop={'size': 6})

    ax2S.plot(x_s_a, z1_s_a, 'b+', markersize = 3, label = 'auto')
    ax2S.plot(x_s_m, z1_s_m, 'r+', markersize = 3, label = 'manual')
    ax2S.set_ylabel('PA (degree) with superficial aponeurosis', fontsize= 8)
    ax2S.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    ax2S.legend(loc = 'upper left', prop={'size': 6})

    ax3S.plot(x_s_a, z2_s_a, 'b+', markersize = 3, label = 'auto')
    ax3S.plot(x_s_m, z2_s_m, 'r+', markersize = 3, label = 'manual')
    ax3S.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    ax3S.set_xlabel('Distance from upper left corner of image (mm)', fontsize= 8)
    ax3S.legend(loc = 'upper left', prop={'size': 6})    


    #panoramic
    figP.suptitle('Panoramic images outputs for participant: ' + participant)

    ax1P.plot(x_p_a, y_p_a, 'b+', markersize = 3, label = 'auto')
    ax1P.plot(x_p_m, y_p_m, 'r+', markersize = 3, label = 'manual')
    ax1P.set_ylabel('FL (mm)', fontsize= 8)
    ax1P.set_xlabel('Distance from insertion (mm)', fontsize= 8)
    ax1P.legend(loc = 'upper left', prop={'size': 6})

    ax2P.plot(x_p_a, z1_p_a, 'b+', markersize = 3, label = 'auto')
    ax2P.plot(x_p_m, z1_p_m, 'r+', markersize = 3, label = 'manual')
    ax2P.set_ylabel('PA (degree) with superficial apoenurosis', fontsize= 8)
    ax2P.set_xlabel('Distance from insertion (mm)', fontsize= 8)
    ax2P.legend(loc = 'upper left', prop={'size': 6})

    ax3P.plot(x_p_a, z2_p_a, 'b+', markersize = 3, label = 'auto')
    ax3P.plot(x_p_m, z2_p_m, 'r+', markersize = 3, label = 'manual')
    ax3P.set_ylabel('PA (degree) with deep aponeurosis', fontsize= 8)
    ax3P.set_xlabel('Distance from insertion (mm)', fontsize= 8)
    ax3P.legend(loc = 'upper left', prop={'size': 6})
    
    plt.show()

    
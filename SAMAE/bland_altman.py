#execute plotBA function to plot all Bland-Altman plots (one per parameter per format of image)
#figures are saved in the folder where path2savefolder leads

import matplotlib.pyplot as plt

def generateType1(title, titleX, titleY):
    fig, ax = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25))
    fig.suptitle(title, fontsize=40, va = 'bottom')
    ax.set_ylabel(titleY, fontsize= 55)
    ax.set_xlabel(titleX, fontsize= 55)
    ax.grid(True)
    ax.tick_params(axis='x', labelsize=35, colors='k')
    ax.tick_params(axis='y', labelsize=35, colors='k')
    return fig, ax


def blandaltman(figureAxis, listManual, listAuto, UoM):
    """Create a Bland-Altman plot from an
    existing figure axis figureAxis, a list with
    the first analysis method (listManual) and a 
    list with the second analysismethod (listAuto)
    
    Reference:
    Chapter 204 - Bland-Altman Plot Analysis, NCSS.com, pp. 204-7: 204-9
    """
    import numpy as np
    import math as m
    
    ####
    #SUBJECT DIFFERENCES METHOD
    listL = np.array([len(l) for l in listAuto])
    listD = np.array([np.array(listAuto[i])-np.array(listManual[i]) for i in range(len(listL))])
    meanD_Subject = np.array([np.mean(l) for l in listD])
    meanD = np.sum(meanD_Subject)/len(meanD_Subject)
    varD_Subject = np.array([np.sum((listD[j]-meanD_Subject[j])**2/(listL[j]-1)) for j in range(len(listL))])
    
    #Within subject random error - this is a squared standard deviation
    Sdw2 = np.sum(varD_Subject*(listL-1)/(np.sum(listL)-len(listL)))
    #Between subject random error - this is a squared standard deviation
    Sdb2 = np.sum((meanD_Subject-meanD)**2)/(len(listL)-1)
    #Harmonic mean of the replicate counts
    mh=len(listL)/np.sum(1/listL)
    #Standard deviation of a difference - this is a squared standard deviation
    Sd2 = Sdb2 + (1-1/mh)*Sdw2
    #Limits of agreement
    LoA_lower = meanD - 1.96 * m.sqrt(Sd2)
    LoA_upper = meanD + 1.96 * m.sqrt(Sd2)
    #95% Confidence interval for LoA - delta method
    v = Sdb2/len(listL) + 1.96**2 / (2* Sd2) * (Sdb2**2/(len(listL)-1) + (1- 1/mh)**2 * Sdw2**2/(np.sum(listL)-len(listL)))
    intLoA_upper = [LoA_upper - 1.96 * m.sqrt(v), LoA_upper + 1.96 * m.sqrt(v)]
    intLoA_lower = [LoA_lower - 1.96 * m.sqrt(v), LoA_lower + 1.96 * m.sqrt(v)]
    
    # VISUALIZATION
    listM = [item for sublist in listManual for item in sublist]
    listA = [item for sublist in listAuto for item in sublist]
    mean = (np.array(listA) + np.array(listM))/2
    diff = np.array(listA) - np.array(listM)
    figureAxis.plot(mean, diff, color = 'k', marker='.', markersize = 40, linestyle='None')
    #mean and limits of agreement axes:
    figureAxis.axhline(meanD, color=(1,0,0), linestyle='--', linewidth=10, label = 'Mean: '+str(round(meanD,2))+UoM)
    figureAxis.axhline(LoA_lower, color='gray', linestyle=':', linewidth=10, label = 'Lower LoA: '+str(round(LoA_lower,2))+UoM)
    figureAxis.axhline(LoA_upper, color='gray', linestyle=':', linewidth=10, label = 'Upper LoA: '+str(round(LoA_upper,2))+UoM)
    #confidence intervals of limits of agreement:
    (lim1, lim2) = figureAxis.get_xbound()
    figureAxis.axhspan(intLoA_upper[0], intLoA_upper[1], color = (0,1,0,0.2))
    figureAxis.axhspan(intLoA_lower[0], intLoA_lower[1], color = (0,1,0,0.2))
    figureAxis.axhspan(intLoA_lower[1], intLoA_upper[0], color = 'gray', alpha=0.2)

    #proportionnal bias via linear regression
    xlim = figureAxis.get_xlim()
    from sklearn.linear_model import LinearRegression
    LinModel = LinearRegression()
    reg = LinModel.fit(np.array(mean).reshape(-1, 1), np.array(diff).reshape(-1, 1))
    a = round(reg.coef_[0][0],2)
    b = round(reg.intercept_[0],2)
    r_sq = round(LinModel.score(np.array(mean).reshape(-1, 1), np.array(diff).reshape(-1, 1)),3)
    Xnew = np.arange(xlim[0], max(mean)+5,1)
    Ynew = LinModel.predict(Xnew.reshape(-1, 1))
    figureAxis.plot(Xnew, Ynew, color = 'k', linestyle = '-', linewidth = 8, label = 'y =' + str(a) + "x + " + str(b)+','+r'$R^{2} =$' + str(r_sq))
    
    box = figureAxis.get_position()
    figureAxis.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
    # Put a legend below current axis   
    figureAxis.legend(loc = 'lower center', bbox_to_anchor = (0.5,-0.3), ncol = 2, fontsize = 35, numpoints=2)
    ###END



def plotBA(path_to_dict, name_dict, path2savefolder=None):
    
    
    participants=['01_Kevin', '02_rafaelopes', '03_charlesbarrand', '04_guilhem',\
        '05_leandre', '06_thomasmartine', '10_victor',\
        '11_youssouf', '12_sufyan', '16_julien', '34_nicolas']
    
    import numpy as np
    import fnmatch
    from dictmanager import load_obj
    
    #NOTES TO UNDERSTAND NOTATIONS
    #d = list of distances from reference point
    #fl = list of FL
    #PAs = list of PA sup
    #PAi = list of PA inf
    #mt = list of MT
    #_s : simple images
    #_p : panoramic images
    #_m : manual
    #_a : automated
    #_filtered: matched fascicles only
    
    participants=['01_Kevin', '02_rafaelopes', '03_charlesbarrand', '04_guilhem',\
        '05_leandre', '06_thomasmartine', '10_victor',\
        '11_youssouf', '12_sufyan', '16_julien', '34_nicolas']

    '************************************************************************'
    '*****************************INITIALIZATION*****************************'
    
    d_s_m = [[] for par in range(len(participants))]
    mt_s_m = [[] for par in range(len(participants))]
    d_s_m_filtered = [[] for par in range(len(participants))]
    fl_s_m_filtered = [[] for par in range(len(participants))]
    PAs_s_m_filtered = [[] for par in range(len(participants))]
    PAi_s_m_filtered = [[] for par in range(len(participants))]

    d_s_a = [[] for par in range(len(participants))]
    mt_s_a = [[] for par in range(len(participants))]
    d_s_a_filtered = [[] for par in range(len(participants))]
    fl_s_a_filtered = [[] for par in range(len(participants))]
    PAs_s_a_filtered = [[] for par in range(len(participants))]
    PAi_s_a_filtered = [[] for par in range(len(participants))]
 
    d_p_m = [[] for par in range(len(participants))]
    mt_p_m = [[] for par in range(len(participants))]
    d_p_m_filtered = [[] for par in range(len(participants))]
    fl_p_m_filtered = [[] for par in range(len(participants))]
    PAs_p_m_filtered = [[] for par in range(len(participants))]
    PAi_p_m_filtered = [[] for par in range(len(participants))]
    
    d_p_a = [[] for par in range(len(participants))]
    mt_p_a = [[] for par in range(len(participants))]
    d_p_a_filtered = [[] for par in range(len(participants))]
    fl_p_a_filtered = [[] for par in range(len(participants))]
    PAs_p_a_filtered = [[] for par in range(len(participants))]
    PAi_p_a_filtered = [[] for par in range(len(participants))]

    #stats on the number of fascicles detected
    nb_fasc_tot_s = 0
    nb_fasc_in_s = 0
    nb_fasc_filt_s = 0
    nb_images_s = 0
    nb_fasc_tot_p = 0
    nb_fasc_in_p = 0
    nb_fasc_filt_p = 0
    nb_images_p = 0

    '************************************************************************'
    '*****************************DATA RETRIEVAL*****************************'
    
    dictio = load_obj(name_dict, path_to_dict)
    l2 = ['fasc*', 'fsc_*']
    
    for par in range(len(participants)):
        
        participant = participants[par]
        fam_folders = [str(d) for d in dictio[participant].keys()]

        s_manuFasc = []
        s_autoFasc = []
        p_manuFasc = []
        p_autoFasc = []
        
        for fam in fam_folders:
            
            ###################################################################
            # simple images
            dictioS = dictio[participant][fam]['BF']['simple']
            images = [str(im) for im in dictioS.keys()]
            for i in images:
                nb_images_s = nb_images_s + 1
                
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
                ###############################################################
                # SIMPLE - automatic
                if ('architecture auto' in dictioS[i]):
                    dictioA = dictioS[i]['architecture auto']
                    midRow = np.mean(dictioA['crop']['lines'])
                    midCol = np.mean(dictioA['crop']['columns'])
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        nb_fasc_tot_s = nb_fasc_tot_s + len(fascicles)
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
                                    nb_fasc_in_s = nb_fasc_in_s + 1
                                
                        if ('MT for labelled points' in dictioM['MT']):
                            for ind0 in range(len(dictioM['MT']['MT for labelled points'])):
                                elem = dictioM['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_s_m[par].append(elem) #MT in mm
                            
                            for ind0 in range(len(dictioA['MT']['MT for labelled points'])):
                                elem = dictioA['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_s_a[par].append(elem)
                                            
            ###################################################################
            # panoramic images
            dictioP = dictio[participant][fam]['BF']['panoramic']
            images = [str(im) for im in dictioP.keys()]
            for i in images:
                nb_images_p = nb_images_p + 1
                
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

                ###############################################################
                # PANORAMIC - automatic
                if ('architecture auto' in dictioP[i]):
                    dictioA = dictioP[i]['architecture auto']
                    if dictioA and ('MT' in dictioA):
                        fascicles = [fa for fa in dictioA if any(fnmatch.fnmatch(fa, p) for p in l2)]
                        nb_fasc_tot_p = nb_fasc_tot_p + len(fascicles)
                        for f in fascicles:
                            dictioF = dictioA[f]
                            idf = fam + '/' + i + '/' + f
                            #only keep fascicles that are entirely within the cropped image,
                            #to compare with manually identified fascicles
                            if len(dictioF.keys())>1 and dictioF['FL']['in/out of the image'] == 'in image':
                                nb_fasc_in_p = nb_fasc_in_p + 1
                                p_autoFasc.append((idf, dictioF['dist from insertion in mm']))
                                d_p_a[par].append(dictioF['dist from insertion in mm'])
                                
                        if ('MT for labelled points' in dictioM['MT']):
                            for ind0 in range(len(dictioM['MT']['MT for labelled points'])):
                                elem = dictioM['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_p_m[par].append(elem) #MT in mm
                            
                            for ind0 in range(len(dictioA['MT']['MT for labelled points'])):
                                elem = dictioA['MT']['MT for labelled points'][ind0]
                                if elem != 'error':
                                    mt_p_a[par].append(elem)

        '************************************************************************'
        '********************MATCHING AUTO & MANUAL FASCICLES*******************'
        
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
        for element in uniqueMatching:
            pathA = element[1].split('/')
            pathM = element[0].split('/')
            nb_fasc_filt_s = nb_fasc_filt_s + 1
            d_s_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['dist from (0,0) of RGB image, in mm'])
            fl_s_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['FL']['length in mm'])
            PAs_s_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['PAsup']['value in degree'])
            PAi_s_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['simple'][pathM[1]]['architecture manual'][pathM[2]]['PAinf']['value in degree'])
            d_s_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['dist from (0,0) of RGB image, in mm'])
            fl_s_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['FL']['length in mm'])
            PAs_s_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['PAsup']['value in degree'])
            PAi_s_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['simple'][pathA[1]]['architecture auto'][pathA[2]]['PAinf']['value in degree'])
        
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
        for element in uniqueMatching:
            pathA = element[1].split('/')
            pathM = element[0].split('/')
            nb_fasc_filt_p = nb_fasc_filt_p + 1
            d_p_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['dist from insertion in mm'])
            fl_p_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['FL']['length in mm'])
            PAs_p_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['PAsup']['value in degree'])
            PAi_p_m_filtered[par].append(dictio[participant][pathM[0]]['BF']['panoramic'][pathM[1]]['architecture manual'][pathM[2]]['PAinf']['value in degree'])
            d_p_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['dist from insertion in mm'])
            fl_p_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['FL']['length in mm'])
            PAs_p_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['PAsup']['value in degree'])
            PAi_p_a_filtered[par].append(dictio[participant][pathA[0]]['BF']['panoramic'][pathA[1]]['architecture auto'][pathA[2]]['PAinf']['value in degree'])

    #generate figures
    figS_fl, axS_fl = generateType1('B-A, FL, simple',r'$mean(FL_a, FL_m) (mm)$', r'$FL_a - FL_m (mm)$')
    figS_pas, axS_pas = generateType1('B-A, PAs, simple',r'$mean(PA_{sup,a} , PA_{sup,m}) (degree)$', r'$PA_{sup,a} - PA_{sup,m} (degree)$')
    figS_pai, axS_pai = generateType1('B-A, PAi, simple',r'$mean(PA_{inf,a} , PA_{inf,m}) (degree)$', r'$PA_{inf,a} - PA_{inf,m} (degree)$')
    figS_mt, axS_mt = generateType1('B-A, MT, simple images', r'$mean(MT_{a} , MT_{m}) (mm)$', r'$MT_{a} - MT_{m} (mm)$')
    figP_fl, axP_fl = generateType1('B-A, FL, pano', r'$mean(FL_a, FL_m) (mm)$', r'$FL_a - FL_m (mm)$')
    figP_pas, axP_pas = generateType1('B-A, PAs, pano',r'$mean(PA_{sup,a} , PA_{sup,m}) (degree)$', r'$PA_{sup,a} - PA_{sup,m} (degree)$')
    figP_pai, axP_pai = generateType1('B-A, PAi, pano', r'$mean(PA_{inf,a} , PA_{inf,m}) (degree)$', r'$PA_{inf,a} - PA_{inf,m} (degree)$')
    figP_mt, axP_mt = generateType1('B-A, MT, panoramic images', r'$mean(MT_{a} , MT_{m}) (mm)$', r'$MT_{a} - MT_{m} (mm)$')

    blandaltman(axS_fl, fl_s_m_filtered, fl_s_a_filtered, UoM = 'mm')
    blandaltman(axP_fl, fl_p_m_filtered, fl_p_a_filtered, UoM = 'mm')    
    blandaltman(figureAxis=axS_pas, listManual=PAs_s_m_filtered, listAuto=PAs_s_a_filtered, UoM=r'$°$')
    blandaltman(axP_pas, PAs_p_m_filtered, PAs_p_a_filtered, UoM = '°')    
    blandaltman(axS_pai, PAi_s_m_filtered, PAi_s_a_filtered, UoM = '°')
    blandaltman(axP_pai, PAi_p_m_filtered, PAi_p_a_filtered, UoM = '°')    
    blandaltman(axS_mt, mt_s_m, mt_s_a, UoM = 'mm')
    blandaltman(axP_mt, mt_p_m, mt_p_a, UoM = 'mm') 
    
    plt.show()
    if path2savefolder is not None:
        figS_fl.savefig('path2savefolder/BA_FLs.jpg')
        figS_pas.savefig('path2savefolder/BA_PASs.jpg')
        figS_pai.savefig('path2savefolder/BA_PAIs.jpg')
        figS_mt.savefig('path2savefolder/BA_MTs.jpg')
        figP_fl.savefig('path2savefolder/BA_FLp.jpg')
        figP_pas.savefig('path2savefolder/BA_PASp.jpg')
        figP_pai.savefig('path2savefolder/BA_PAIp.jpg')
        figP_mt.savefig('path2savefolder/BA_MTp.jpg')


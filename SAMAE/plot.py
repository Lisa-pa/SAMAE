"""plots"""
import matplotlib.pyplot as plt

def rand_jitter(arr, sensib = 0.01, lowerLimit=None, upperLimit=None):
    """
    Creation of jittering in a one-D data array arr
    sensibility can be adjust with parameter sensib.
    In case of an array only made of a same value, you can use the lowerlimit and
    upperlimit parameters
    """
    import numpy as np
    if len(arr) != 0:
        if lowerLimit is None or upperLimit is None:
            stdev = sensib * (np.amax(arr) - np.amin(arr))
        else:
            stdev = sensib * (upperLimit-lowerLimit)
        return arr + np.random.randn(len(arr)) * stdev
    else:
        return []

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
    
    
#PLOTS WITH TYPE FUNCTIONS
def plotType1(figureAxis, listManual, listAuto, c):
        figureAxis.plot(listManual, listAuto, alpha = 0.5, color = c, marker='.', markersize = 40, linestyle='None')
        
def plotType2(figureAxis, listManual, listAuto, listAuto2):
        figureAxis.plot(rand_jitter([0]*len(listAuto), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listAuto, alpha = 0.3, color = (0.8,0.2,1), marker='.', markersize = 40, linestyle='None')
        figureAxis.plot(rand_jitter([2]*len(listManual), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listManual, alpha = 0.3, color = (0.8,0.2,1), marker='.', markersize = 40, linestyle='None')
        figureAxis.plot(rand_jitter([1]*len(listAuto2), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listAuto2, alpha = 0.3, color = (0.8,0.2,1), marker='.', markersize = 40, linestyle='None')

def plotType3(titre_participant, titreY, figureAxis, listAuto, listAuto2, listManual, mean1, mean2, mean3, std1, std2, std3,c):
        figureAxis.plot(rand_jitter([0]*len(listAuto), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listAuto, alpha = 0.5, color = c, marker='.', markersize = 7, linestyle='None')
        figureAxis.plot(rand_jitter([1]*len(listAuto2), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listAuto2, alpha = 0.5, color = c, marker='.', markersize = 7, linestyle='None')
        figureAxis.plot(rand_jitter([2]*len(listManual), sensib=0.07, lowerLimit=-0.45, upperLimit=-0.05), listManual, alpha = 0.5, color = c, marker='.', markersize = 7, linestyle='None')
        figureAxis.set_xticks([0,1,2])
        figureAxis.set_xticklabels(['all auto', 'filtered', 'manual'], rotation=25, ha='right', fontsize=8)    
        figureAxis.errorbar([0.1,1.1,2.1], [mean1, mean2, mean3], [std1,std2,std3], fmt='ok', markersize=4, lw=1, capsize = 3, capthick = 1, barsabove=True)
        figureAxis.tick_params(axis='x', labelsize=8, colors='k')
        figureAxis.tick_params(axis='y', labelsize=8, colors='k')
        figureAxis.grid(True)
        figureAxis.set_title(titre_participant, loc='center', color = c, fontsize = 10)

#GENERATE THE FIGURES TEMPLATE WITH FUNCTIONS
def generateType0(title, titleY):
    fig, ax = plt.subplots(1,1,sharex = False, sharey = False, figsize =(4,7), dpi = 800)
    fig.suptitle(title, va = 'bottom', fontsize = 40)
    ax.set_ylabel(titleY, fontsize= 11)
    ax.grid(True)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['automatic', 'manual'], fontsize=11)
    ax.tick_params(axis='x', colors='k')
    ax.tick_params(axis='y', colors='k')
    return fig, ax

def generateType1(title, titleX, titleY):
    fig, ax = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25), dpi = 800)
    fig.suptitle(title, fontsize=40, va = 'bottom')
    ax.set_ylabel(titleY, fontsize= 55)
    ax.set_xlabel(titleX, fontsize= 55)
    ax.grid(True)
    ax.tick_params(axis='x', labelsize=35, colors='k')
    ax.tick_params(axis='y', labelsize=35, colors='k')
    return fig, ax

def generateType2(title, titleY):
    fig, ax = plt.subplots(1,1,sharex = False, sharey = False, figsize =(25,25), dpi = 800)
    fig.suptitle(title, va = 'bottom', fontsize=40)
    ax.grid(True)
    ax.set_ylabel(titleY, fontsize= 55)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(['all auto', 'filtered', 'manual'], rotation=45, ha='right', fontsize=55)    
    ax.tick_params(axis='x', colors='k')
    ax.tick_params(axis='y', labelsize=35, colors='k')
    return fig, ax

def generateType3(titre, Ytitle):
    fig = plt.figure(figsize =(19/2.54, 35/2.54), dpi = 800)
    ax = fig.add_subplot(111)
    fig.suptitle(titre)
    gs = fig.add_gridspec(4,3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    ax7 = fig.add_subplot(gs[2, 0])
    ax8 = fig.add_subplot(gs[2, 1])
    ax9 = fig.add_subplot(gs[2, 2])
    ax10 = fig.add_subplot(gs[3, 0])
    ax11 = fig.add_subplot(gs[3, 1])
    ax12 = fig.add_subplot(gs[3, 2])
    ax12.set_visible(False)
    l = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]
    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    # Set common labels
    ax.set_ylabel(Ytitle)
    
    return fig,l
    
    
#FUNCTION THAT RETRIEVE DATA FROM RECORDED DICTIONNARY
def plotFeatures(path_to_dict, name_dict, colors, participants):
    """
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
    #_p : panoramic images
    #_m : manual
    #_a : automated
    



    '************************************************************************'
    '*****************************INITIALIZATION*****************************'
    
    
    d_s_m = [[] for par in range(len(participants))]
    fl_s_m = [[] for par in range(len(participants))]
    PAs_s_m = [[] for par in range(len(participants))]
    PAi_s_m = [[] for par in range(len(participants))]
    mt_s_m = [[] for par in range(len(participants))]
    d_s_m_filtered = [[] for par in range(len(participants))]
    fl_s_m_filtered = [[] for par in range(len(participants))]
    PAs_s_m_filtered = [[] for par in range(len(participants))]
    PAi_s_m_filtered = [[] for par in range(len(participants))]

    d_s_a = [[] for par in range(len(participants))]
    fl_s_a = [[] for par in range(len(participants))]
    PAs_s_a = [[] for par in range(len(participants))]
    PAi_s_a = [[] for par in range(len(participants))]
    mt_s_a = [[] for par in range(len(participants))]
    d_s_a_filtered = [[] for par in range(len(participants))]
    fl_s_a_filtered = [[] for par in range(len(participants))]
    PAs_s_a_filtered = [[] for par in range(len(participants))]
    PAi_s_a_filtered = [[] for par in range(len(participants))]
 
    d_p_m = [[] for par in range(len(participants))]
    fl_p_m = [[] for par in range(len(participants))]
    PAs_p_m = [[] for par in range(len(participants))]
    PAi_p_m = [[] for par in range(len(participants))]
    mt_p_m = [[] for par in range(len(participants))]
    d_p_m_filtered = [[] for par in range(len(participants))]
    fl_p_m_filtered = [[] for par in range(len(participants))]
    PAs_p_m_filtered = [[] for par in range(len(participants))]
    PAi_p_m_filtered = [[] for par in range(len(participants))]
    
    d_p_a = [[] for par in range(len(participants))]
    fl_p_a = [[] for par in range(len(participants))]
    PAs_p_a = [[] for par in range(len(participants))]
    PAi_p_a = [[] for par in range(len(participants))]
    mt_p_a = [[] for par in range(len(participants))]
    d_p_a_filtered = [[] for par in range(len(participants))]
    fl_p_a_filtered = [[] for par in range(len(participants))]
    PAs_p_a_filtered = [[] for par in range(len(participants))]
    PAi_p_a_filtered = [[] for par in range(len(participants))]

    diff_calfct_p = []
    diff_calfct_s = []

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
                                    fl_s_a[par].append(dictioF['FL']['length in mm'])
                                    PAs_s_a[par].append(dictioF['PAsup']['value in degree'])
                                    PAi_s_a[par].append(dictioF['PAinf']['value in degree'])
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
                                            
                    
                    #calibration factors difference
                    diff_calfct_s.append(abs(dictioM['calfct_to_mm']-\
                                             dictioA['calfct_to_mm']['horizontal axis']))
            

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
                        fl_p_m[par].append(dictioF['FL']['length in mm'])
                        PAs_p_m[par].append(dictioF['PAsup']['value in degree'])
                        PAi_p_m[par].append(dictioF['PAinf']['value in degree'])            

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
    

    '************************************************************************'
    '**************************PLOTS INITIALIZATION**************************'
    
#    #plots for simple images
#    # -- calibration
#    figSc, (axSc1, axSc2) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15), dpi = 800)
#    # -- MT
#    figS0, axS0 = generateType0('Mean MT per participant in simple images', 'MT (mm)')
#    # -- FL
#    figS1, axS1 = generateType1('Automatically estimated length of fascicles relative to their closest manually identified fascicle, in simple images for all participants', r'$FL_m (mm)$', r'$FL_a (mm)$')
#    figS2, axS2 = generateType2('Comparison of automatic and manual outputs for FL estimation, in simple images for all participants', 'FL (mm)')
#    figS3, Slist_axes_FL = generateType3('Comparison of automatic and manual outputs for FL estimation, in simple images per participant', 'FL (mm)')
#    figS4, axS4 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of fascicles length, in simple images',r'$mean(FL_a, FL_m) (mm)$', r'$FL_a - FL_m (mm)$')

#    # -- PAs
#    figS5, axS5 = generateType1('Automatically estimated superior PA of fascicles relative to their closest manually identified fascicle, in simple images for all participants', r'$PA_{sup,m} (degree)$', r'$PA_{sup,a} (degree)$')
#    figS6, axS6 = generateType2('Comparison of automatic and manual outputs for superior PA estimation, in simple images for all participants', r'$PA_{sup} (degree)$')
#    figS7, Slist_axes_PAs = generateType3('Comparison of automatic and manual outputs for superior PA estimation, in simple images per participant', r'$PA_{sup}$ (degree)')
#    figS8, axS8 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of superior PA, in simple images',r'$mean(PA_{sup,a} , PA_{sup,m}) (degree)$', r'$PA_{sup,a} - PA_{sup,m} (degree)$')

#    # -- PAi
#    figS9, axS9 = generateType1('Automatically estimated inferior PA of fascicles relative to their closest manually identified fascicle, in simple images for all participants', r'$PA_{inf,m} (degree)$', r'$PA_{inf,a} (degree)$')
#    figS10, axS10 = generateType2('Comparison of automatic and manual outputs for inferior PA estimation, in simple images for all participants', r'$PA_{inf} (degree)$')
#    figS11, Slist_axes_PAi = generateType3('Comparison of automatic and manual outputs for inferior PA estimation, in simple images per participant', r'$PA_{inf}$ (degree)')
#    figS12, axS12 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of inferior PA, in simple images',r'$mean(PA_{inf,a} , PA_{inf,m}) (degree)$', r'$PA_{inf,a} - PA_{inf,m} (degree)$')
#
#    #plots for panoramic images
#    # -- calibration
#    figPc, (axPc1, axPc2) = plt.subplots(1,2,sharex = False, sharey = False, figsize =(15,15), dpi = 800)

#    # -- MT
#    figP0, axP0 = generateType0('Mean MT per participant in panoramic images','MT (mm)')

#    # -- FL
#    figP1, axP1 =generateType1('Automatically estimated length of fascicles relative to their closest manually identified fascicle, in panoramic images for all participants', r'$FL_m (mm)$', r'$FL_a (mm)$')
#    figP2, axP2 = generateType2('Comparison of automatic and manual outputs for FL estimation, in panoramic images for all participants', 'FL (mm)')   
#    figP3, Plist_axes_FL = generateType3('Comparison of automatic and manual outputs for FL estimation, in panoramic images per participant','FL (mm)')   
#    figP4, axP4 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of fascicles length, in panoramic images', r'$mean(FL_a, FL_m) (mm)$', r'$FL_a - FL_m (mm)$')

#    # -- PAs
#    figP5, axP5 =generateType1('Automatically estimated superior PA of fascicles relative to their closest manually identified fascicle, in panoramic images for all participants', r'$PA_{sup,m} (degree)$', r'$PA_{sup,a} (degree)$')
#    figP6, axP6 = generateType2('Comparison of automatic and manual outputs for superior PA estimation, in panoramic images for all participants', r'$PA_{sup} (degree)$')
#    figP7, Plist_axes_PAs = generateType3('Comparison of automatic and manual outputs for superior PA estimation, in panoramic images per participant', r'$PA_{sup}$ (degree)')
#    figP8, axP8 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of superior PA, in panoramic images',r'$mean(PA_{sup,a} , PA_{sup,m}) (degree)$', r'$PA_{sup,a} - PA_{sup,m} (degree)$')

#    # -- PAi
#    figP9, axP9 = generateType1('Automatically estimated inferior PA of fascicles relative to their closest manually identified fascicle, in panoramic images for all participants', r'$PA_{inf,m} (degree)$', r'$PA_{inf,a} (degree)$')
#    figP10, axP10 = generateType2('Comparison of automatic and manual outputs for inferior PA estimation, in panoramic images for all participants', r'$PA_{inf} (degree)$')
#    
#    figP11, Plist_axes_PAi = generateType3('Comparison of automatic and manual outputs for inferior PA estimation, in panoramic images per participant', r'$PA_{inf}$ (degree)')
#    
#    figP12, axP12 = generateType1('Bland-Altman plot for the comparison of manual and automatic estimation of inferior PA, in panoramic images', r'$mean(PA_{inf,a} , PA_{inf,m}) (degree)$', r'$PA_{inf,a} - PA_{inf,m} (degree)$')
    
    
    for par in range(len(participants)):
        ID = participants[par][:2]


        '**********************************************************************'
        '******************************STATISTICS******************************'

        #statistics per participant
        
        p_a_PASUP_std = np.std(PAs_p_a[par])
        p_a_PASUP_mean = np.mean(PAs_p_a[par])
        p_a_filt_PASUP_std = np.std(PAs_p_a_filtered[par])
        p_a_filt_PASUP_mean = np.mean(PAs_p_a_filtered[par])
        p_a_PAINF_std = np.std(PAi_p_a[par])
        p_a_PAINF_mean = np.mean(PAi_p_a[par]) 
        p_a_filt_PAINF_std = np.std(PAi_p_a_filtered[par])
        p_a_filt_PAINF_mean = np.mean(PAi_p_a_filtered[par]) 
        p_a_FL_std = np.std(fl_p_a[par])
        p_a_FL_mean = np.mean(fl_p_a[par])
        p_a_filt_FL_std = np.std(fl_p_a_filtered[par])
        p_a_filt_FL_mean = np.mean(fl_p_a_filtered[par])
        p_a_MT_mean = np.mean(mt_p_a[par])
        p_a_MT_std =np.std(mt_p_a[par])
        #*****#
        p_m_PASUP_std = np.std(PAs_p_m[par])
        p_m_PASUP_mean = np.mean(PAs_p_m[par])
        p_m_PAINF_std = np.std(PAi_p_m[par])
        p_m_PAINF_mean = np.mean(PAi_p_m[par])
        p_m_FL_std = np.std(fl_p_m[par])
        p_m_FL_mean = np.mean(fl_p_m[par])
        p_m_MT_mean = np.mean(mt_p_m[par])
        p_m_MT_std = np.std(mt_p_m[par])
        p_m_filt_PASUP_std =np.std(PAs_p_m_filtered[par])
        p_m_filt_PASUP_mean =np.mean(PAs_p_m_filtered[par])
        p_m_filt_PAINF_std =np.std(PAi_p_m_filtered[par])
        p_m_filt_PAINF_mean =np.mean(PAi_p_m_filtered[par])        
        p_m_filt_FL_std =np.std(fl_p_m_filtered[par])
        p_m_filt_FL_mean =np.mean(fl_p_m_filtered[par])
        #*****#
        s_a_PASUP_std =np.std(PAs_s_a[par])
        s_a_PASUP_mean =np.mean(PAs_s_a[par])
        s_a_PAINF_std =np.std(PAi_s_a[par])
        s_a_PAINF_mean =np.mean(PAi_s_a[par])        
        s_a_FL_std =np.std(fl_s_a[par])
        s_a_FL_mean =np.mean(fl_s_a[par])
        s_a_MT_mean = np.mean(mt_s_a[par])
        s_a_MT_std =np.std(mt_s_a[par])
        s_a_filt_PASUP_std =np.std(PAs_s_a_filtered[par])
        s_a_filt_PASUP_mean =np.mean(PAs_s_a_filtered[par])
        s_a_filt_PAINF_std =np.std(PAi_s_a_filtered[par])
        s_a_filt_PAINF_mean =np.mean(PAi_s_a_filtered[par])        
        s_a_filt_FL_std =np.std(fl_s_a_filtered[par])
        s_a_filt_FL_mean =np.mean(fl_s_a_filtered[par])
        #*****#
        s_m_PASUP_std =np.std(PAs_s_m[par])
        s_m_PASUP_mean =np.mean(PAs_s_m[par])
        s_m_PAINF_std =np.std(PAi_s_m[par])
        s_m_PAINF_mean =np.mean(PAi_s_m[par])
        s_m_FL_std =np.std(fl_s_m[par])
        s_m_FL_mean =np.mean(fl_s_m[par])
        s_m_MT_mean = np.mean(mt_s_m[par])
        s_m_MT_std = np.std(mt_s_m[par])
        s_m_filt_PASUP_std =np.std(PAs_s_m_filtered[par])
        s_m_filt_PASUP_mean =np.mean(PAs_s_m_filtered[par])
        s_m_filt_PAINF_std =np.std(PAi_s_m_filtered[par])
        s_m_filt_PAINF_mean =np.mean(PAi_s_m_filtered[par])        
        s_m_filt_FL_std =np.std(fl_s_m_filtered[par])
        s_m_filt_FL_mean =np.mean(fl_s_m_filtered[par])
        #*****#
        
        """print('STATISTICS PER PARTICIPANT')
        print('part = ', ID, 'stats=',\
              p_a_PASUP_std, p_a_PASUP_mean,\
              p_a_PAINF_std, p_a_PAINF_mean, p_a_FL_std, p_a_FL_mean,\
              p_a_filt_PASUP_std, p_a_filt_PASUP_mean,\
              p_a_filt_PAINF_std, p_a_filt_PAINF_mean, p_a_filt_FL_std, p_a_filt_FL_mean,\
              p_a_MT_mean, p_a_MT_std,\
              p_m_PASUP_std, p_m_PASUP_mean,\
              p_m_PAINF_std, p_m_PAINF_mean, p_m_FL_std, p_m_FL_mean,\
              p_m_filt_PASUP_std, p_m_filt_PASUP_mean,\
              p_m_filt_PAINF_std, p_m_filt_PAINF_mean, p_m_filt_FL_std, p_m_filt_FL_mean,\
              p_m_MT_mean, p_m_MT_std,\
              s_a_PASUP_std, s_a_PASUP_mean,\
              s_a_PAINF_std, s_a_PAINF_mean, s_a_FL_std, s_a_FL_mean,\
              s_a_filt_PASUP_std, s_a_filt_PASUP_mean,\
              s_a_filt_PAINF_std, s_a_filt_PAINF_mean, s_a_filt_FL_std, s_a_filt_FL_mean,\
              s_a_MT_mean, s_a_MT_std,\
              s_m_PASUP_std, s_m_PASUP_mean,\
              s_m_PAINF_std, s_m_PAINF_mean, s_m_FL_std, s_m_FL_mean,\
              s_m_filt_PASUP_std, s_m_filt_PASUP_mean,\
              s_m_filt_PAINF_std, s_m_filt_PAINF_mean, s_m_filt_FL_std, s_m_filt_FL_mean,\
              s_m_MT_mean, s_m_MT_std)"""
                
        '*********************************************************************'
        '********************************PLOTS********************************'

#        ##########
#        # -- FL
#        # --- simple
#        # ---- auto = f(manu) with filtered fascicles of all participants
#        plotType1(axS1, fl_s_m_filtered[par], fl_s_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT; with jitter and alpha value
#        plotType2(axS2, fl_s_m_filtered[par], fl_s_a[par], fl_s_a_filtered[par])
#        # ---- PER PARTICIPANT; with jitter and alpha value
#        plotType3(ID, 'FL (mm)', Slist_axes_FL[par], fl_s_a[par], fl_s_a_filtered[par], fl_s_m_filtered[par], s_a_FL_mean, s_a_filt_FL_mean, s_m_filt_FL_mean, s_a_FL_std, s_a_filt_FL_std, s_m_filt_FL_std,colors[par])
#                
#        # --- panoramic
#        # ---- auto = f(manu)
#        plotType1(axP1, fl_p_m_filtered[par], fl_p_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT
#        plotType2(axP2, fl_p_m_filtered[par], fl_p_a[par], fl_p_a_filtered[par])
#        # ---- PER PARTICIPANT
#        plotType3(ID, 'FL (mm)', Plist_axes_FL[par], fl_p_a[par], fl_p_a_filtered[par], fl_p_m_filtered[par], p_a_FL_mean, p_a_filt_FL_mean, p_m_filt_FL_mean, p_a_FL_std, p_a_filt_FL_std, p_m_filt_FL_std,colors[par])
#
#        ##########
#        # -- PA sup
#        # --- simple
#        # ---- auto = f(manu)
#        plotType1(axS5, PAs_s_m_filtered[par], PAs_s_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT; with jitter and alpha value
#        plotType2(axS6, PAs_s_m_filtered[par], PAs_s_a[par], PAs_s_a_filtered[par])
#        # ---- PER PARTICIPANT; with jitter and alpha value
#        plotType3(ID, r'$PA_{sup} (degree)$', Slist_axes_PAs[par], PAs_s_a[par], PAs_s_a_filtered[par], PAs_s_m_filtered[par], s_a_PASUP_mean, s_a_filt_PASUP_mean, s_m_filt_PASUP_mean, s_a_PASUP_std, s_a_filt_PASUP_std, s_m_filt_PASUP_std, colors[par])
#        
#        # --- panoramic
#        # ---- auto = f(manu)
#        plotType1(axP5, PAs_p_m_filtered[par], PAs_p_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT
#        plotType2(axP6, PAs_p_m_filtered[par], PAs_p_a[par], PAs_p_a_filtered[par])
#        # ---- PER PARTICIPANT
#        plotType3(ID, r'$PA_{sup} (degree)', Plist_axes_PAs[par], PAs_p_a[par], PAs_p_a_filtered[par], PAs_p_m_filtered[par], p_a_PASUP_mean, p_a_filt_PASUP_mean, p_m_filt_PASUP_mean, p_a_PASUP_std, p_a_filt_PASUP_std, p_m_filt_PASUP_std, colors[par])
#        
#        ##########
#        # -- PA inf
#        # --- simple
#        # ---- auto = f(manu)
#        plotType1(axS9, PAi_s_m_filtered[par], PAi_s_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT; with jitter and alpha value
#        plotType2(axS10, PAi_s_m_filtered[par], PAi_s_a[par], PAi_s_a_filtered[par])
#        # ---- PER PARTICIPANT; with jitter and alpha value
#        plotType3(ID, r'$PA_{inf} (degree)', Slist_axes_PAi[par], PAi_s_a[par], PAi_s_a_filtered[par], PAi_s_m_filtered[par], s_a_PAINF_mean, s_a_filt_PAINF_mean, s_m_filt_PAINF_mean, s_a_PAINF_std, s_a_filt_PAINF_std, s_m_filt_PAINF_std, colors[par])
#        
#        # --- panoramic
#        # ---- auto = f(manu)
#        plotType1(axP9, PAi_p_m_filtered[par], PAi_p_a_filtered[par], colors[par])
#        # ---- TOT PARTICIPANT
#        plotType2(axP10, PAi_p_m_filtered[par], PAi_p_a[par], PAi_p_a_filtered[par])
#        # ---- PER PARTICIPANT
#        plotType3(ID, r'$PA_{inf} (degree)', Plist_axes_PAi[par], PAi_p_a[par], PAi_p_a_filtered[par], PAi_p_m_filtered[par], p_a_PAINF_mean, p_a_filt_PAINF_mean, p_m_filt_PAINF_mean, p_a_PAINF_std, p_a_filt_PAINF_std, p_m_filt_PAINF_std, colors[par])
#                
#        ##########  
#
#        #MT 
#        axS0.plot(rand_jitter([0,1]), rand_jitter([s_a_MT_mean, s_m_MT_mean]), color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 1)
#        axS0.text(0.3+par*0.3/11, (s_a_MT_mean+s_m_MT_mean)/2-0.2, ID, color = colors[par], fontsize = 8)
#        axP0.plot(rand_jitter([0,1]), rand_jitter([p_a_MT_mean, p_m_MT_mean]), color = colors[par], marker='.', markersize = 7, linestyle='solid', linewidth = 1)
#        axP0.text(0.3+par*0.3/11, (p_a_MT_mean+p_m_MT_mean)/2-0.2, ID, color = colors[par], fontsize = 8)
    
    #statistics on the population
    TOT_p_a_PASUP_mean = np.mean([item for sublist in PAs_p_a for item in sublist])
    TOT_p_a_PAINF_mean = np.mean([item for sublist in PAi_p_a for item in sublist])
    TOT_p_a_FL_mean = np.mean([item for sublist in fl_p_a for item in sublist])
    TOT_p_a_filt_PASUP_mean = np.mean([item for sublist in PAs_p_a_filtered for item in sublist])
    TOT_p_a_filt_PAINF_mean = np.mean([item for sublist in PAi_p_a_filtered for item in sublist])
    TOT_p_a_filt_FL_mean = np.mean([item for sublist in fl_p_a_filtered for item in sublist])
    TOT_p_a_MT_mean = np.mean([item for sublist in mt_p_a for item in sublist])
    TOT_p_a_PASUP_std = np.std([item for sublist in PAs_p_a for item in sublist])
    TOT_p_a_PAINF_std = np.std([item for sublist in PAi_p_a for item in sublist])
    TOT_p_a_FL_std = np.std([item for sublist in fl_p_a for item in sublist])
    TOT_p_a_filt_PASUP_std = np.std([item for sublist in PAs_p_a_filtered for item in sublist])
    TOT_p_a_filt_PAINF_std = np.std([item for sublist in PAi_p_a_filtered for item in sublist])
    TOT_p_a_filt_FL_std = np.std([item for sublist in fl_p_a_filtered for item in sublist])
    TOT_p_a_MT_std =np.std([item for sublist in mt_p_a for item in sublist])
    
    TOT_p_m_PASUP_mean = np.mean([item for sublist in PAs_p_m for item in sublist])
    TOT_p_m_PAINF_mean = np.mean([item for sublist in PAi_p_m for item in sublist])
    TOT_p_m_FL_mean = np.mean([item for sublist in fl_p_m for item in sublist])
    TOT_p_m_MT_mean = np.mean([item for sublist in mt_p_m for item in sublist])
    TOT_p_m_PASUP_std = np.std([item for sublist in PAs_p_m for item in sublist])
    TOT_p_m_PAINF_std = np.std([item for sublist in PAi_p_m for item in sublist])
    TOT_p_m_FL_std = np.std([item for sublist in fl_p_m for item in sublist])
    TOT_p_m_MT_std = np.std([item for sublist in mt_p_m for item in sublist])
    TOT_p_m_filt_PASUP_mean =np.mean([item for sublist in PAs_p_m_filtered for item in sublist])
    TOT_p_m_filt_PAINF_mean =np.mean([item for sublist in PAi_p_m_filtered for item in sublist])
    TOT_p_m_filt_FL_mean =np.mean([item for sublist in fl_p_m_filtered for item in sublist])
    TOT_p_m_filt_PASUP_std =np.std([item for sublist in PAs_p_m_filtered for item in sublist])
    TOT_p_m_filt_PAINF_std =np.std([item for sublist in PAi_p_m_filtered for item in sublist])
    TOT_p_m_filt_FL_std =np.std([item for sublist in fl_p_m_filtered for item in sublist])

    TOT_s_a_PASUP_mean =np.mean([item for sublist in PAs_s_a for item in sublist])
    TOT_s_a_PAINF_mean =np.mean([item for sublist in PAi_s_a for item in sublist])
    TOT_s_a_FL_mean =np.mean([item for sublist in fl_s_a for item in sublist])
    TOT_s_a_filt_PASUP_mean =np.mean([item for sublist in PAs_s_a_filtered for item in sublist])
    TOT_s_a_filt_PAINF_mean =np.mean([item for sublist in PAi_s_a_filtered for item in sublist])
    TOT_s_a_filt_FL_mean =np.mean([item for sublist in fl_s_a_filtered for item in sublist])
    TOT_s_a_MT_mean = np.mean([item for sublist in mt_s_a for item in sublist])
    TOT_s_a_PASUP_std =np.std([item for sublist in PAs_s_a for item in sublist])
    TOT_s_a_PAINF_std =np.std([item for sublist in PAi_s_a for item in sublist])
    TOT_s_a_FL_std =np.std([item for sublist in fl_s_a for item in sublist])
    TOT_s_a_filt_PASUP_std =np.std([item for sublist in PAs_s_a_filtered for item in sublist])
    TOT_s_a_filt_PAINF_std =np.std([item for sublist in PAi_s_a_filtered for item in sublist])
    TOT_s_a_filt_FL_std =np.std([item for sublist in fl_s_a_filtered for item in sublist])
    TOT_s_a_MT_std =np.std([item for sublist in mt_s_a for item in sublist])

    TOT_s_m_PASUP_mean =np.mean([item for sublist in PAs_s_m for item in sublist])
    TOT_s_m_PAINF_mean =np.mean([item for sublist in PAi_s_m for item in sublist])
    TOT_s_m_FL_mean =np.mean([item for sublist in fl_s_m for item in sublist])
    TOT_s_m_MT_mean = np.mean([item for sublist in mt_s_m for item in sublist])
    TOT_s_m_PASUP_std =np.std([item for sublist in PAs_s_m for item in sublist])
    TOT_s_m_PAINF_std =np.std([item for sublist in PAi_s_m for item in sublist])
    TOT_s_m_FL_std =np.std([item for sublist in fl_s_m for item in sublist])
    TOT_s_m_MT_std = np.std([item for sublist in mt_s_m for item in sublist])
    TOT_s_m_filt_PASUP_mean =np.mean([item for sublist in PAs_s_m_filtered for item in sublist])
    TOT_s_m_filt_PAINF_mean =np.mean([item for sublist in PAi_s_m_filtered for item in sublist])
    TOT_s_m_filt_FL_mean =np.mean([item for sublist in fl_s_m_filtered for item in sublist])
    TOT_s_m_filt_PASUP_std =np.std([item for sublist in PAs_s_m_filtered for item in sublist])
    TOT_s_m_filt_PAINF_std =np.std([item for sublist in PAi_s_m_filtered for item in sublist])
    TOT_s_m_filt_FL_std =np.std([item for sublist in fl_s_m_filtered for item in sublist])
    
    print('ALL PARTICIPANTS STATISTICS')
    print('pano auto PAS: ', TOT_p_a_filt_PASUP_mean, '+/-', TOT_p_a_filt_PASUP_std,\
    'pano auto PAi: ', TOT_p_a_filt_PAINF_mean, '+/-', TOT_p_a_filt_PAINF_std,\
    'pano auto FL: ', TOT_p_a_filt_FL_mean, '+/-', TOT_p_a_filt_FL_std,\
    'pano manu PAs: ',TOT_p_m_filt_PASUP_mean, '+/-', TOT_p_m_filt_PASUP_std,\
    'pano manu PAi: ', TOT_p_m_filt_PAINF_mean, TOT_p_m_filt_PAINF_std,\
    'pano manu FL: ', TOT_p_m_filt_FL_mean, TOT_p_m_filt_FL_std,\
    'simple auto PAS: ', TOT_s_a_filt_PASUP_mean, '+/-', TOT_s_a_filt_PASUP_std,\
    'simple auto PAi: ', TOT_s_a_filt_PAINF_mean, '+/-', TOT_s_a_filt_PAINF_std,\
    'simple auto FL: ', TOT_s_a_filt_FL_mean, '+/-', TOT_s_a_filt_FL_std,\
    'simple manu PAs: ',TOT_s_m_filt_PASUP_mean, '+/-', TOT_s_m_filt_PASUP_std,\
    'simple manu PAi: ', TOT_s_m_filt_PAINF_mean, TOT_s_m_filt_PAINF_std,\
    'simple manu FL: ', TOT_s_m_filt_FL_mean, TOT_s_m_filt_FL_std)

    #Plot error bars
        #on MT plots
#    axP0.errorbar([0,1],[TOT_p_a_MT_mean, TOT_p_m_MT_mean], [TOT_p_a_MT_std,TOT_p_m_MT_std], fmt='ok', lw=2, capsize = 5, capthick = 2)
#    axS0.errorbar([0,1],[TOT_s_a_MT_mean, TOT_s_m_MT_mean], [TOT_s_a_MT_std,TOT_s_m_MT_std], fmt='ok', lw=2, capsize = 5, capthick = 2)
#        # on FL plot for all participants
#    axS2.errorbar([0.15,1.15,2.15], [TOT_s_a_FL_mean, TOT_s_a_filt_FL_mean, TOT_s_m_filt_FL_mean], [TOT_s_a_FL_std, TOT_s_a_filt_FL_std, TOT_s_m_filt_FL_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)
#    axP2.errorbar([0.15,1.15,2.15], [TOT_p_a_FL_mean, TOT_p_a_filt_FL_mean, TOT_p_m_filt_FL_mean], [TOT_p_a_FL_std, TOT_p_a_filt_FL_std, TOT_p_m_filt_FL_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)
#        # on PAs plot for all participants
#    axS6.errorbar([0.15,1.15,2.15], [TOT_s_a_PASUP_mean, TOT_s_a_filt_PASUP_mean, TOT_s_m_filt_PASUP_mean], [TOT_s_a_PASUP_std, TOT_s_a_filt_PASUP_std, TOT_s_m_filt_PASUP_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)
#    axP6.errorbar([0.15,1.15,2.15], [TOT_p_a_PASUP_mean, TOT_p_a_filt_PASUP_mean, TOT_p_m_filt_PASUP_mean], [TOT_p_a_PASUP_std, TOT_p_a_filt_PASUP_std, TOT_p_m_filt_PASUP_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)
#        # on PAi plot for all participants
#    axS10.errorbar([0.15,1.15,2.15], [TOT_s_a_PAINF_mean, TOT_s_a_filt_PAINF_mean, TOT_s_m_filt_PAINF_mean], [TOT_s_a_PAINF_std, TOT_s_a_filt_PAINF_std, TOT_s_m_filt_PAINF_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)
#    axP10.errorbar([0.15,1.15,2.15], [TOT_p_a_PAINF_mean, TOT_p_a_filt_PAINF_mean, TOT_p_m_filt_PAINF_mean], [TOT_p_a_PAINF_std, TOT_p_a_filt_PAINF_std, TOT_p_m_filt_PAINF_std], fmt='ok', markersize = 30, lw=8, capsize = 16, capthick = 8, barsabove=True)

    #Bland-Altman plots
    # - FL
        # ---- Bland-Altman
#    blandaltman(axS4, fl_s_m_filtered, fl_s_a_filtered, UoM = 'mm')
#    blandaltman(axP4, fl_p_m_filtered, fl_p_a_filtered, UoM = 'mm')    
#    blandaltman(figureAxis=axS8, listManual=PAs_s_m_filtered, listAuto=PAs_s_a_filtered, UoM=r'$°$')
#    blandaltman(axP8, PAs_p_m_filtered, PAs_p_a_filtered, UoM = '°')    
#    blandaltman(axS12, PAi_s_m_filtered, PAi_s_a_filtered, UoM = '°')
#    blandaltman(axP12, PAi_p_m_filtered, PAi_p_a_filtered, UoM = '°')    
   


#    #calibration factors
#    if len(diff_calfct_s)>0:
#        diff_med_s = np.median(diff_calfct_s)
#        diff_max_s = max(diff_calfct_s)
#        diff_mean_s = np.mean(diff_calfct_s)
#        diff_std_s = np.std(diff_calfct_s)
#        figSc.suptitle('Difference in manual/automatic calibration factor in simple images', va = 'bottom')
#        axSc1.plot([1]*len(diff_calfct_s), diff_calfct_s, color = colors[1], marker = 'o', markersize = 5, linestyle = 'None')
#        axSc1.set_ylabel(r'$ |Calib_m - Calib_a| (mm/pixel)$', fontsize= 8)
#        axSc1.grid(True)
#        axSc2.text(0.2,0.2, 'Mean ='+str(diff_mean_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axSc2.text(0.2,0.7, 'Median ='+str(diff_med_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axSc2.text(0.2,0.5, 'Max ='+str(diff_max_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axSc2.text(0.2,0.1, 'STD ='+str(diff_std_s)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axSc2.set_axis_off()
#
#    if len(diff_calfct_p) >0:
#        diff_med_p = np.median(diff_calfct_p)
#        diff_max_p = max(diff_calfct_p)
#        diff_mean_p = np.mean(diff_calfct_p)
#        diff_std_p = np.std(diff_calfct_p)
#        figPc.suptitle('Difference in manual/automatic calibration factor in panoramic images',  va = 'bottom')
#        axPc1.plot([1]*len(diff_calfct_p), diff_calfct_p, color = colors[1], marker = 'o', markersize = 5, linestyle = 'None')
#        axPc1.set_ylabel(r'$ |Calib_m - Calib_a| (mm/pixel)$', fontsize= 8)
#        axPc1.grid(True)
#        axPc2.text(0.2,0.2, 'Mean ='+str(diff_mean_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axPc2.text(0.2,0.7, 'Median ='+str(diff_med_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axPc2.text(0.2,0.5, 'Max ='+str(diff_max_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axPc2.text(0.2,0.1, 'STD ='+str(diff_std_p)+' mm/pixel', fontsize= 15, color='k', bbox=dict(facecolor='red', alpha=0.5))
#        axPc2.set_axis_off()

#    figS13, axS13 = generateType1('B-A plot, MT, simple images', r'$mean(MT_{a} , MT_{m}) (mm)$', r'$MT_{a} - MT_{m} (mm)$')
#    figP13, axP13 = generateType1('B-A plot, MT, panoramic images', r'$mean(MT_{a} , MT_{m}) (mm)$', r'$MT_{a} - MT_{m} (mm)$')
#    blandaltman(axS13, mt_s_m, mt_s_a, UoM = 'mm')
#    blandaltman(axP13, mt_p_m, mt_p_a, UoM = 'mm') 
#    figS13.savefig('C:/Users/Lisa Paillard/Desktop/BA_MTs.jpg')
#    figP13.savefig('C:/Users/Lisa Paillard/Desktop/BA_MTp.jpg')

#    figSc.savefig('C:/Users/Lisa Paillard/Desktop/calibs.jpg')
#    figS0.savefig('C:/Users/Lisa Paillard/Desktop/mts.jpg')
#    figS1.savefig('C:/Users/Lisa Paillard/Desktop/fls1.jpg')
#    figS2.savefig('C:/Users/Lisa Paillard/Desktop/fls2.jpg')
#    figS3.savefig('C:/Users/Lisa Paillard/Desktop/fls3.jpg')
#    figS4.savefig('C:/Users/Lisa Paillard/Desktop/fls4.jpg')
#    figS5.savefig('C:/Users/Lisa Paillard/Desktop/pasups1.jpg')
#    figS6.savefig('C:/Users/Lisa Paillard/Desktop/pasups2.jpg')
#    figS7.savefig('C:/Users/Lisa Paillard/Desktop/pasups3.jpg')
#    figS8.savefig('C:/Users/Lisa Paillard/Desktop/pasups4.jpg')
#    figS9.savefig('C:/Users/Lisa Paillard/Desktop/painfs1.jpg')
#    figS10.savefig('C:/Users/Lisa Paillard/Desktop/painfs2.jpg')
#    figS11.savefig('C:/Users/Lisa Paillard/Desktop/painfs3.jpg')
#    figS12.savefig('C:/Users/Lisa Paillard/Desktop/painfs4.jpg')
#    figPc.savefig('C:/Users/Lisa Paillard/Desktop/calibp.jpg')
#    figP0.savefig('C:/Users/Lisa Paillard/Desktop/mtp.jpg')
#    figP1.savefig('C:/Users/Lisa Paillard/Desktop/flp1.jpg')
#    figP2.savefig('C:/Users/Lisa Paillard/Desktop/flp2.jpg')
#    figP3.savefig('C:/Users/Lisa Paillard/Desktop/flp3.jpg')
#    figP4.savefig('C:/Users/Lisa Paillard/Desktop/flp4.jpg')
#    figP5.savefig('C:/Users/Lisa Paillard/Desktop/pasupp1.jpg')
#    figP6.savefig('C:/Users/Lisa Paillard/Desktop/pasupp2.jpg')
#    figP7.savefig('C:/Users/Lisa Paillard/Desktop/pasupp3.jpg')
#    figP8.savefig('C:/Users/Lisa Paillard/Desktop/pasupp4.jpg')
#    figP9.savefig('C:/Users/Lisa Paillard/Desktop/painfp1.jpg')
#    figP10.savefig('C:/Users/Lisa Paillard/Desktop/painfp2.jpg')
#    figP11.savefig('C:/Users/Lisa Paillard/Desktop/painfp3.jpg')
#    figP12.savefig('C:/Users/Lisa Paillard/Desktop/painfp4.jpg')
    
    #stats on detection of fascicles
    #average number of detected fascicles per image
    meanF_tot_s = nb_fasc_tot_s / nb_images_s
    meanF_tot_p = nb_fasc_tot_p / nb_images_p
    #average number of fascicles IN per image
    meanF_in_s = nb_fasc_in_s / nb_images_s
    meanF_in_p = nb_fasc_in_p / nb_images_p
    #average number for the comparison with manual data
    meanF_filtered_s = nb_fasc_filt_s / nb_images_s
    meanF_filtered_p = nb_fasc_filt_p / nb_images_p
    
    print('Statistics on fascicles detection')
    print('Average nb of fasc detected in SI:', meanF_tot_s)
    print('Average nb of fasc detected in PI:', meanF_tot_p)
    print('Average nb of IN fasc detected in SI:', meanF_in_s)
    print('Average nb of IN fasc detected in PI:', meanF_in_p)
    print('Average nb of filtered fasc detected in SI:', meanF_filtered_s)
    print('Average nb of filtered fasc detected in PI:', meanF_filtered_p)
    print('Total number of filtered fasc detected in SI:',nb_fasc_filt_s)
    print('Total number of filtered fasc detected in PI:',nb_fasc_filt_p)
    
    nb_filt_fascPerParts = [len(subL) for subL in fl_s_a_filtered]
    nb_filt_fascPerPartp = [len(subL) for subL in fl_p_a_filtered]
    print('s', np.mean(nb_filt_fascPerParts))
    print('p', np.mean(nb_filt_fascPerPartp))
#    from scipy.stats.mstats import ttest_rel
#    t,p=ttest_rel([item for sublist in PAs_s_m_filtered for item in sublist],[item for sublist in PAs_s_a_filtered for item in sublist],axis=None)
#    print('PAS s',t,p)
#    t2,p2=ttest_rel([item for sublist in PAs_p_m_filtered for item in sublist],[item for sublist in PAs_p_a_filtered for item in sublist],axis=None)
#    print('PAS p',t2,p2)
#    t3,p3=ttest_rel([item for sublist in PAi_s_m_filtered for item in sublist],[item for sublist in PAi_s_a_filtered for item in sublist],axis=None)
#    print('PAI s',t3,p3)
#    t4,p4=ttest_rel([item for sublist in PAi_p_m_filtered for item in sublist],[item for sublist in PAi_p_a_filtered for item in sublist],axis=None)
#    print('PAI p',t4,p4)
#    t5,p5=ttest_rel([item for sublist in fl_s_m_filtered for item in sublist],[item for sublist in fl_s_a_filtered for item in sublist],axis=None)
#    print('FL s',t5,p5)
#    t6,p6=ttest_rel([item for sublist in fl_p_m_filtered for item in sublist],[item for sublist in fl_p_a_filtered for item in sublist],axis=None)
#    print('FL p',t6,p6)
#    t7,p7=ttest_rel([item for sublist in mt_s_m for item in sublist],[item for sublist in mt_s_a for item in sublist],axis=None)
#    print('mt s',t7,p7)
#    t8,p8=ttest_rel([item for sublist in mt_p_m for item in sublist],[item for sublist in mt_p_a for item in sublist],axis=None)
#    print('mt p',t8,p8)
    
    ###
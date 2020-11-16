def sizeEffect(listManual, listAuto):
    import numpy as np
    import math as m
    listM = [item for sublist in listManual for item in sublist]
    listA = [item for sublist in listAuto for item in sublist]
    diff = np.array(listA) - np.array(listM)
    meanD = np.mean(diff)
    stdD = np.std(diff)
    #Cohen size effect
    size_effect = meanD/stdD
    cohen_d = (np.mean(listA)-np.mean(listM)) / (m.sqrt((np.std(listA) ** 2 + np.std(listM) ** 2) / 2))
    return size_effect, cohen_d
    
def ttests(path_to_dict, name_dict='\\TOTresults'):
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
                if par == 9 and fam =='fam_2' and i=='img_2':
                    print(par, fam, i)
                else:
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

    
    #t_tests
    print('paired samples t-tests resuts: ')
    from scipy.stats.mstats import ttest_rel
    #NOTES: we cannot user '..._filtered' arrays directly because of their structure
    #we need to flatten them to 1-D lists
    t,p=ttest_rel([item for sublist in PAs_s_m_filtered for item in sublist],[item for sublist in PAs_s_a_filtered for item in sublist],axis=None)
    print('PAS s',p)
    t2,p2=ttest_rel([item for sublist in PAs_p_m_filtered for item in sublist],[item for sublist in PAs_p_a_filtered for item in sublist],axis=None)
    print('PAS p',p2)
    t3,p3=ttest_rel([item for sublist in PAi_s_m_filtered for item in sublist],[item for sublist in PAi_s_a_filtered for item in sublist],axis=None)
    print('PAI s',p3)
    t4,p4=ttest_rel([item for sublist in PAi_p_m_filtered for item in sublist],[item for sublist in PAi_p_a_filtered for item in sublist],axis=None)
    print('PAI p',p4)
    t5,p5=ttest_rel([item for sublist in fl_s_m_filtered for item in sublist],[item for sublist in fl_s_a_filtered for item in sublist],axis=None)
    print('FL s',p5)
    t6,p6=ttest_rel([item for sublist in fl_p_m_filtered for item in sublist],[item for sublist in fl_p_a_filtered for item in sublist],axis=None)
    print('FL p',p6)
    t7,p7=ttest_rel([item for sublist in mt_s_m for item in sublist],[item for sublist in mt_s_a for item in sublist],axis=None)
    print('mt s',p7)
    t7_2,p7_2=ttest_rel([np.mean(sublist) for sublist in mt_s_m],[np.mean(sublist) for sublist in mt_s_a],axis=None)
    print('mt s for means',p7_2)
    t8,p8=ttest_rel([item for sublist in mt_p_m for item in sublist],[item for sublist in mt_p_a for item in sublist],axis=None)
    print('mt p',p8)
    t8_2,p8_2=ttest_rel([np.mean(sublist) for sublist in mt_p_m],[np.mean(sublist) for sublist in mt_p_a],axis=None)
    print('mt p for means',p8_2)

    print('independent samples t-tests resuts: ')
    from scipy.stats.mstats import ttest_rel, ttest_ind
    #NOTES: we cannot user '..._filtered' arrays directly because of their structure
    #we need to flatten them to 1-D lists
    t,p=ttest_ind([item for sublist in PAs_s_m_filtered for item in sublist],[item for sublist in PAs_s_a_filtered for item in sublist],axis=None)
    print('PAS s',p)
    t2,p2=ttest_ind([item for sublist in PAs_p_m_filtered for item in sublist],[item for sublist in PAs_p_a_filtered for item in sublist],axis=None)
    print('PAS p',p2)
    t3,p3=ttest_ind([item for sublist in PAi_s_m_filtered for item in sublist],[item for sublist in PAi_s_a_filtered for item in sublist],axis=None)
    print('PAI s',p3)
    t4,p4=ttest_ind([item for sublist in PAi_p_m_filtered for item in sublist],[item for sublist in PAi_p_a_filtered for item in sublist],axis=None)
    print('PAI p',p4)
    t5,p5=ttest_ind([item for sublist in fl_s_m_filtered for item in sublist],[item for sublist in fl_s_a_filtered for item in sublist],axis=None)
    print('FL s',p5)
    t6,p6=ttest_ind([item for sublist in fl_p_m_filtered for item in sublist],[item for sublist in fl_p_a_filtered for item in sublist],axis=None)
    print('FL p',p6)
    t7,p7=ttest_ind([item for sublist in mt_s_m for item in sublist],[item for sublist in mt_s_a for item in sublist],axis=None)
    print('mt s',p7)
    t7_2,p7_2=ttest_ind([np.mean(sublist) for sublist in mt_s_m],[np.mean(sublist) for sublist in mt_s_a],axis=None)
    print('mt s for means',p7_2)
    t8,p8=ttest_ind([item for sublist in mt_p_m for item in sublist],[item for sublist in mt_p_a for item in sublist],axis=None)
    print('mt p',p8)
    t8_2,p8_2=ttest_ind([np.mean(sublist) for sublist in mt_p_m],[np.mean(sublist) for sublist in mt_p_a],axis=None)
    print('mt p for means',p8_2)


    
    #size effects
    s1=sizeEffect(PAs_s_m_filtered,PAs_s_a_filtered)
    s2=sizeEffect(PAs_p_m_filtered,PAs_p_a_filtered)
    s3=sizeEffect(PAi_s_m_filtered,PAi_s_a_filtered)
    s4=sizeEffect(PAi_p_m_filtered,PAi_p_a_filtered)
    s5=sizeEffect(fl_s_m_filtered,fl_s_a_filtered)
    s6=sizeEffect(fl_p_m_filtered,fl_p_a_filtered)
    s7=sizeEffect(mt_s_m,mt_s_a)
    s8=sizeEffect(mt_p_m,mt_p_a)
    print('Size effects: ')
    print('PAS s',s1)
    print('PAS p',s2)
    print('PAi s',s3)
    print('PAi p',s4)
    print('fl s',s5)
    print('fl p',s6)
    print('mt s',s7)
    print('mt p',s8)

    mt_s_a_filt = [[] for par in range(len(participants))]
    mt_s_m_filt = [[] for par in range(len(participants))]
    for p in range(len(mt_s_a)):
        for val in range(len(mt_s_a[p])):
            if p==9:
                if mt_s_a[p][val]>mt_s_m[p][val]+2.37 or mt_s_a[p][val]<mt_s_m[p][val]-2.08:
                    print('aberrante valeur: participant ', p, ' , place val ', val)
                else:
                    mt_s_a_filt[p].append(mt_s_a[p][val])
                    mt_s_m_filt[p].append(mt_s_m[p][val])
            else:
                mt_s_a_filt[p].append(mt_s_a[p][val])
                mt_s_m_filt[p].append(mt_s_m[p][val])

    print('apres avoir enleve les valeurs out of LoA: ')
    t7,p7=ttest_rel([item for sublist in mt_s_m_filt for item in sublist],[item for sublist in mt_s_a_filt for item in sublist],axis=None)
    print('mt s',p7)
        
        
        
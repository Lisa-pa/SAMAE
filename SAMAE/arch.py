"""Management of the dict architecture"""

def dame_arch_paths(path_to_folders, participants, test='architecture'):
    """Create dict of paths were architecture data will be stored
    """

    import os
    import fnmatch

    arch_paths = dict()
    for participante in range(len(participants)):
        l1 = os.listdir(str(path_to_folders) + '\\' + participants[participante])
        l2 = ['POST', '*fam_*']
        fam_folders = [x for x in l1 if any(fnmatch.fnmatch(x, p) for p in l2)]

        # create list fo inside fam folders
        arch_paths[participants[participante]] = {}
        for fa in range(len(fam_folders)):
            ll = os.listdir(
                str(path_to_folders) + '\\' + str(participants[participante]) + '\\' + str(fam_folders[fa]))
            if test in ll:
                arch_paths[participants[participante]][fam_folders[fa]] = str(
                    path_to_folders) + '\\' + str(participants[participante]) + '\\' + str(fam_folders[fa]) + '\\' + str(test)

    return arch_paths


def dame_arch_data(archpaths):
    """""This function gathers architecture data from .txt files
    
    Arguments:
        archpaths {list} -- list of paths to every folder containing architecture data
        
    Returns:
        [dict] -- Returns dictionary containing all data for each path
    """

    import os
    import fnmatch

    archdata = dict()
    for participant in archpaths.keys():
        
        sessions = dict()
        for ntest in archpaths[participant].keys():
            print('Este es el test de: ', participant, 'sesion: ', ntest)
            path = archpaths[participant][ntest]

            # '_bfs.txt' and '_sms.txt' refer to single image files only
            files = os.listdir(path)

            bf_files = {'simple': fnmatch.filter(files, str('*_bfs.txt')),
                        'panoramic': fnmatch.filter(files, str('*_bfp.txt')),
                        'landmark': fnmatch.filter(files, str('*_bfpm.txt'))}
            '''
            sm_files = {'simple': fnmatch.filter(files, str('*_sms.txt')),
                        'panoramic': fnmatch.filter(files, str('*_smp.txt')),
                        'landmark': fnmatch.filter(files, str('*_smpm.txt'))}
            '''
            # get data
            bf_data = _organize_arch(fils=bf_files, pth=path)
            '''
            sm_data = _organize_arch(fils=sm_files, pth=path)
            '''
            
            muscles = dict()
            muscles['BF'] = bf_data
            '''
            muscles['SM'] = sm_data
            '''
            
            sessions[ntest] = muscles
            
        archdata[participant] = sessions
    

    return _calcula_arch(data=archdata)


def _calcula_arch(data):
    """Compute manual and automatic architectural features from 
    dict of coordinates. It updates the input dict
    
    Arguments:
        data {dict} -- dictionary containing coordinates for each subject/trial/image
    
    Returns:
        dict -- Dict containing coordinates and architecture results for each participant/trial/image
    """

    import numpy as np
    import autoP as autoP
    import autoS as autoS
    import manuP as manuP
    import manuS as manuS
    
    for part in data.keys():
        for ntest in data[part].keys():
            for msc in data[part][ntest].keys():
                for echo in data[part][ntest][msc].keys():

                    if echo == 'landmark':
                        distances = []
                        for img in data[part][ntest][msc][echo].keys():
                            a = data[part][ntest][msc]['landmark'][img]['coords']
                            distances.append(distsimpleimg(coords=a))
                        data[part][ntest][msc]['landmark']['dist'] = np.mean(distances)
                    
                    
                    for img in data[part][ntest][msc][echo].keys():

                        print(part, ntest, msc, echo, img)

                        if echo == 'panoramic':
                            path_to_txt = data[part][ntest][msc][echo][img]['path']
                            path_to_jpg = path_to_txt[:-3] + 'jpg'                            
                            idx = data[part][ntest][msc]['panoramic'][img]['coords']
                            # manual processing
                            architecture1 = idfascicles(coords=idx, img='pano')
                            points_sup_m = architecture1['aposup']['coords']
                            points_inf_m = architecture1['apoinf']['coords']
                            architecture1 = manuP.panoManu(architecture1, None)
                            data[part][ntest][msc]['panoramic'][img]['architecture manual'] = architecture1
                            #automatic processing
                            architecture2 = autoP.panoprocessing(path_to_jpg, path_to_txt)
                            if architecture2:
                                points_sup_a = architecture2['aposup']['coords']
                                points_inf_a = architecture2['apoinf']['coords']
                                data[part][ntest][msc]['panoramic'][img]['architecture auto'] = architecture2
                                
                                if ('MT' in architecture2):
                                    MT1, MT2 = forMTcomparison(points_sup_m, points_inf_m, points_sup_a, points_inf_a,architecture1['calfct_to_mm'],\
                                                        architecture2['calfct_to_mm before resize']['vertical axis'])
                                    data[part][ntest][msc]['panoramic'][img]['architecture auto']['MT']['MT for labelled points'] = MT2
                                    data[part][ntest][msc]['panoramic'][img]['architecture manual']['MT']['MT for labelled points'] = MT1

                        if echo == 'simple':
                            path_to_txt = data[part][ntest][msc][echo][img]['path']
                            path_to_jpg = path_to_txt[:-3] + 'jpg'                            
                            idx = data[part][ntest][msc]['simple'][img]['coords']
                            #manual processing
                            architecture1 = idfascicles(coords=idx, img='simple')
                            points_sup_m = architecture1['aposup']['coords']
                            points_inf_m = architecture1['apoinf']['coords']
                            architecture1 = manuS.simpleManu(architecture1)
                            data[part][ntest][msc]['simple'][img]['architecture manual'] = architecture1
                            #automatic processing
                            architecture2 = autoS.simpleprocessing(path_to_jpg)
                            if architecture2:
                                points_sup_a = architecture2['aposup']['coords']
                                points_inf_a = architecture2['apoinf'] ['coords']                           
                                data[part][ntest][msc]['simple'][img]['architecture auto'] = architecture2
                            
                                if ('MT' in architecture2):
                                    MT1, MT2 = forMTcomparison(points_sup_m, points_inf_m, points_sup_a, points_inf_a,\
                                                        architecture1['calfct_to_mm'], architecture2['calfct_to_mm']['vertical axis'])
                                    data[part][ntest][msc]['simple'][img]['architecture auto']['MT']['MT for labelled points'] = MT2
                                    data[part][ntest][msc]['simple'][img]['architecture manual']['MT']['MT for labelled points'] = MT1
    return data


def _organize_arch(fils, pth):

    """Allocate data from each specific type of file (keys from the input dict) to a new dict
    
    Arguments:
        fils {dict} -- Dictionary containing type of files and list of files
    
    Returns:
        [dict] -- [description]
    """

    import numpy as np
    
    imgdata = dict()
    for i in fils.keys():
        
        images = dict()
        for ii in np.arange(len(fils[i])):
            images[str('img_' + str(ii+1))] = {'path': pth + str('\\') + str(fils[i][ii]),
                                               'coords': np.loadtxt(pth + str('\\') + str(fils[i][ii]), skiprows=1, usecols=(-2, -1))}

        imgdata[i] = images    

    return imgdata



def idfascicles(coords, img):
    """Identify coordinates of each fascicle manually labelled in the panoramic (10 points for each fascicle) 
    and the simple images (4 points for each)
    
    Arguments:
        coords {array} -- Array containing `x` and `y` coordinates of points defining scale, 
        insertion point, aponeuroses and fascicles
        
        img {str} -- string (`pano` or `simple`) indicating wether it is a panoramic 
        or simple ultrasound image scan 
    
    Returns:
        dict() -- Dictionary containing the corresponding coordinates for each architectural feature
    """

    import numpy as np

    architecture = dict()
    architecture['calfct_to_mm'] = 10./(coords[1, 1] - coords[0, 1])  # 10 mm calibration factor for all panoramic and simple images

    if img == 'pano':
        coordS = coords[8:13, :]
        coordI = coords[3:8, :]
        # inversion of coordinates, so that points appear as [row, column] and not [column, row]
        coordS[:,[0,1]] = coordS[:,[1,0]]
        coordI[:,[0,1]] = coordI[:,[1,0]]
        #first coordinate is row, second coordinate is column
        architecture['insertion'] = {'coords': [coords[2,1],coords[2,0]]}
        architecture['aposup'] = {'coords': coordS}
        architecture['apoinf'] = {'coords': coordI}

        fsc_idx = np.arange(start=13, stop=len(coords), step=10)
        for f in np.arange(len(fsc_idx)):
            coordF = coords[fsc_idx[f]:fsc_idx[f]+10, :]
            coordF[:,[0,1]] = coordF[:,[1,0]]
            #first coordinate is row, second coordinate is column
            architecture['fsc_'+str(f+1)] = {'coords': coordF}

    if img == 'simple':
        coordS = coords[2:6, :]
        coordS[:,[0,1]] = coordS[:,[1,0]]
        coordI = coords[6:10, :]
        coordI[:,[0,1]] = coordI[:,[1,0]]
        #first coordinate is row, second coordinate is column
        architecture['aposup'] = {'coords': coordS}
        architecture['apoinf'] = {'coords': coordI}

        fsc_idx = np.arange(start=10, stop=len(coords), step=4)
        for f in np.arange(len(fsc_idx)):
            coordF = coords[fsc_idx[f]:fsc_idx[f]+4, :]
            coordF[:,[0,1]] = coordF[:,[1,0]]
            #first coordinate is row, second coordinate is column
            architecture['fsc_'+str(f+1)] = {'coords': coordF}

    return architecture


def distsimpleimg(coords):
    """Return distance (mm) from muscle insertion point 
    to the place where the single image ultrasound was taken
    
    Arguments:
        coords {array} -- Array containing x and y coordinates of calibration scale, 
        insertion point and place of the image.
    
    Returns:
        [float] -- distance in cm
    """

    cal_fct = coords[1, 1] - coords[0, 1]
    insertion = coords[2]
    img_place = coords[3]

    distance = [(insertion - img_place) / cal_fct][0][0]

    return distance

def forMTcomparison(ptsSup_m, ptsInf_m, ptsSup_a, ptsInf_a, calibV_m, calibV_a):
    """
    Computes MT at discrete points to compare labelled data and automatic processing.
    Inputs
        ptsSup_m: list of points of superficial aponeurosis, manual labelling
        ptsInf_m: list of points of deep aponeurosis manual labelling
        ptsSup_a: list of continuous points of superficial aponeurosis, automatic processing, all along the US image
        ptsInf_a: list of continuous points of deep aponeurosis, automatic processing, all along the US image
        calibV_m (float): vertical calibration factor computed with manual labelling
        calibV_a (float): vertical calibration factor computed from automatic processing
    Outputs:
        2 lists of same length, that contains MT from manual data, and MT from 
            automatic processing at the same columns as for manual MT
    """
    MT_m = []
    MT_a = []
    
    for ind in range(len(ptsSup_m)):
        # extract column and rows of manual points
        col_1 = ptsSup_m[ind][1]
        col_2 = ptsInf_m[ind][1]
        lig_1 = ptsSup_m[ind][0]
        lig_2 = ptsInf_m[ind][0]
        lig_3 = -1
        lig_4 = -1
        # compute manual MT
        MT_m.append(abs(lig_2 - lig_1)*calibV_m)
            
        if ptsSup_a != 'error' and ptsInf_a != 'error':
            # look for automatic points that have the same column as manual points
            for ind2 in range(len(ptsSup_a)):
                if int(ptsSup_a[ind2][1]) == int(col_1):
                    lig_3 = ptsSup_a[ind2][0]
                if int(ptsInf_a[ind2][1]) == int(col_2):
                    lig_4 = ptsInf_a[ind2][0]
            
            #compute automatic MT if the previous automatic points were found
            if lig_3>= 0 and lig_4 >=0:
                MT_a.append(abs(lig_4 - lig_3)*calibV_a)
            else:
                MT_a.append('error')
                
    return MT_m, MT_a
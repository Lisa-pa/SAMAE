def dame_participants():
    """Create list containing folder names

    Returns:
        list: list of strings
    """    
    part = ['01_Kevin','02_rafaelopes']
    
    '''
    part = ['01_Kevin', '02_rafaelopes', '03_charlesbarrand', '04_guilhem',
        '05_leandre', '06_thomasmartine', '10_victor', 
        '11_youssouf', '12_sufyan', '16_julien', '34_nicolas']
    '''
    
    # part = ['01_Kevin', '02_rafaelopes', '03_charlesbarrand', '04_guilhem',
    #         '05_leandre', '06_thomasmartine', '09_serge', '10_victor',
    #         '11_youssouf', '12_sufyan', '14_thomasFrancois',
    #         '15_davidsimeon', '16_julien',

    #         '18_mehdi', '20_tiavina', '22_granska', '23_jeandruais',
    #         '24_guillaumeleroy', '25_mouaad', '25_mouaad', '26_waid',
    #         '28_hugo', '29_johann', '30_nelsonlopez', '31_romain',
    #         '41_pierrelouisgauche', '33_gabrielbeq', '34_nicolas',

    #         '07_jihad', '08_jason', '17_simonavr', '19_enzo', '21_vincent',
    #         '27_thibaultrigal', '35_benjamin', '36_baptiste', '37_samuel',
    #         '38_samir', '40_jamesjeremy', '42_javi']

    return part


def dame_arch_paths(path_to_folders, test='architecture'):
    """Create dict of paths were architecture data have been stored

    """

    import os
    import fnmatch

    participants = dame_participants()

    arch_paths = dict()
    for participante in range(len(participants)):
        l1 = os.listdir(str(path_to_folders) + '\\' + participants[participante])
        l2 = ['POST', '*fam_*']
        fam_folders = [x for x in l1 if any(fnmatch.fnmatch(x, p) for p in l2)]

        # create list fo inside fam folders
        arch_paths[participants[participante]] = {}
        for fa in range(len(fam_folders)):
            ll = os.listdir(
                str(path_to_folders) + '\\' + participants[participante] + '\\' + str(fam_folders[fa]))
            if test in ll:
                arch_paths[participants[participante]][fam_folders[fa]] = str(
                    path_to_folders) + '\\' + participants[participante] + '\\' + fam_folders[fa] + '\\' + test

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
    """COmpute architectural features from dict of coordinates. It updates the input dict
    
    Arguments:
        data {dict} -- dictionary containing coordinates for each subject/trial/image
    
    Returns:
        dict -- Dict containing coordinates and architecture results for each participant/trial/image
    """

    import numpy as np
    import autoP
    import autoS
    import manuP
    import manuS
    
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
                            architecture1 = manuP.panoManu(architecture1, None)
                            data[part][ntest][msc]['panoramic'][img]['architecture manual'] = architecture1
                            #automatic processing
                            architecture2 = autoP.panoprocessing(path_to_jpg, path_to_txt)
                            data[part][ntest][msc]['panoramic'][img]['architecture auto'] = architecture2
                            


                        if echo == 'simple':
                            path_to_txt = data[part][ntest][msc][echo][img]['path']
                            path_to_jpg = path_to_txt[:-3] + 'jpg'                            
                            idx = data[part][ntest][msc]['simple'][img]['coords']
                            #manual processing
                            architecture1 = idfascicles(coords=idx, img='simple')
                            architecture1 = manuS.simpleManu(architecture1)
                            data[part][ntest][msc]['simple'][img]['architecture manual'] = architecture1
                            #automatic processing
                            architecture2 = autoS.simpleprocessing(path_to_jpg)
                            data[part][ntest][msc]['simple'][img]['architecture auto'] = architecture2
                            
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
        architecture['insertion'] = {'coords': coords[2,:]}
        architecture['aposup'] = {'coords': coords[8:13, :]}
        architecture['apoinf'] = {'coords': coords[3:8, :]}

        fsc_idx = np.arange(start=13, stop=len(coords), step=10)
        for f in np.arange(len(fsc_idx)):
            architecture['fsc_'+str(f+1)] = {'coords': coords[fsc_idx[f]:fsc_idx[f]+10, :]}

    if img == 'simple':
        architecture['aposup'] = {'coords': coords[2:6, :]}
        architecture['apoinf'] = {'coords': coords[6:10, :]}

        fsc_idx = np.arange(start=10, stop=len(coords), step=4)
        for f in np.arange(len(fsc_idx)):
            architecture['fsc_'+str(f+1)] = {'coords': coords[fsc_idx[f]:fsc_idx[f]+4, :]}

    return architecture

def distsimpleimg(coords):
    """Return distance (cm) from muscle insertion point 
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
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:53:09 2020

@author: Lisa Paillard
"""

def panoManu(manu_pts_dict, path_img = None):
    """
    """
    import MUFeaM
    import numpy as np
    import FaDe
    import fnmatch
    import scipy.interpolate as interpolate

    calib = manu_pts_dict['calfct_to_mm']
    pt_inter = manu_pts_dict['insertion']['coords']

    #aponeuroses interpolation with Bspline
    list_sup = list(manu_pts_dict['aposup']['coords'])
    list_sup.sort(key=lambda x:x[1])
    y_sup = [list_sup[i][1] for i in range(len(list_sup))]
    x_sup = [list_sup[i][0] for i in range(len(list_sup))]
    
    list_inf = list(manu_pts_dict['apoinf']['coords'])
    list_inf.sort(key=lambda x:x[1])
    y_inf = [list_inf[i][1] for i in range(len(list_inf))]
    x_inf = [list_inf[i][0] for i in range(len(list_inf))]
    
    spline_sup = interpolate.UnivariateSpline(y_sup, x_sup, k=1, ext = 0)
    spline_inf = interpolate.UnivariateSpline(y_inf, x_inf, k=1, ext = 0)
    coordSup, spline_sup = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_sup, interval = [0, int(pt_inter[1])])
    coordInf, spline_inf = MUFeaM.pointsCoordinates(typeA = 'spline', param = spline_inf, interval = [0, int(pt_inter[1])])


    # features computation
    abscissa, thickness, splineT = MUFeaM.muscleThickness(start = max(min(y_sup), min(y_inf)), end = min(max(y_sup), max(y_inf)), calibV = calib, calibH = calib, spl1 = spline_sup, spl2 = spline_inf)
    manu_pts_dict['MT'] = {'coords': thickness, 'columns interval': [max(min(y_sup), min(y_inf)), min(max(y_sup), max(y_inf))]}
    
    for fsc in fnmatch.filter(manu_pts_dict.keys(), 'fsc_*'):
        fasci = manu_pts_dict[fsc]['coords']        
        
        #Bspline
        spl = FaDe.approximateFasc(typeapprox = 'Bspline', listF = [fasci], d = 1)
        #find intersection points with aponeuroses
        inters_inf, inters_sup, spl = MUFeaM.findIntersections(spl_inf = spline_inf, spl_sup = spline_sup,\
                                                               listSpl = spl, search_interval = [0, pt_inter[1]], start=0)
        if len(spl)>0:
            #find pennation angles
            PA_sup = MUFeaM.pennationAngles(spl_a = spline_sup, listS_f = spl, listI = inters_sup, xcalib = calib, ycalib =  calib)
            PA_inf = MUFeaM.pennationAngles(spl_a = spline_inf, listS_f = spl, listI = inters_inf, xcalib = calib, ycalib = calib)
            #calculate fascicle length
            FL = MUFeaM.fasciclesLength(listS_f = spl, listIu = inters_sup, listId = inters_inf, xcalib = calib, ycalib = calib)
            loc = MUFeaM.locateFasc(inters_inf, refPoint = pt_inter, ycalib = calib)
            
            manu_pts_dict[fsc]['dist from insertion in mm'] = loc[0]
            manu_pts_dict[fsc]['PAsup'] = {'value in degree' : PA_sup[0],
                         'intersection with apo': inters_sup[0]}
            manu_pts_dict[fsc]['PAinf'] = {'value in degree' : PA_inf[0],
                         'intersection with apo': inters_inf[0]}
            manu_pts_dict[fsc]['FL'] = {'length in mm': FL[0]}

        if path_img is not None: 
            import cv2
            #fascicles visualization
            I = cv2.imread(path_img, -1)
            newy = np.arange(0, I.shape[1], 1)
            newx = np.int32(spl[0](newy))
            coord = np.vstack((newx, newy)).T
            for n2 in range(coord.shape[0]):
                if coord[n2][0]>=0 and coord[n2][0]<I.shape[0]:
                    I[coord[n2][0], coord[n2][1], :] = [0,255,0]
                    if coord[n2][1] - 1 >= 0 and coord[n2][1] - 1 < I.shape[1]:
                        I[coord[n2][0], coord[n2][1] - 1, :] = [0,255,0]
                    if coord[n2][1] + 1 >= 0 and coord[n2][1] + 1 < I.shape[1]:
                        I[coord[n2][0], coord[n2][1] + 1, :] = [0,255,0]

    if path_img is not None:
        # aponeuroses visualization
        for index in range(coordSup.shape[0]):
            if coordSup[index][0] >= 0 and coordSup[index][0] < I.shape[0]:
                I[coordSup[index][0], coordSup[index][1], :] = [255, 0, 0]
                if coordSup[index][0] +1 >= 0 and coordSup[index][0]+1 < I.shape[0]:
                    I[coordSup[index][0]+1, coordSup[index][1],:] = [255, 0, 0]
                if coordSup[index][0]-1 >= 0 and coordSup[index][0]-1 < I.shape[0]:
                    I[coordSup[index][0]-1, coordSup[index][1],:] = [255, 0, 0]
            if coordInf[index][0] >= 0 and coordInf[index][0] < I.shape[0]:
                I[coordInf[index][0], coordInf[index][1], :] = [255, 0, 0]
                if coordInf[index][0]+1 >= 0 and coordInf[index][0]+1 < I.shape[0]:
                    I[coordInf[index][0]+1, coordInf[index][1], :] = [255, 0, 0]
                if coordInf[index][0]-1 >= 0 and coordInf[index][0]-1 < I.shape[0]:
                    I[coordInf[index][0]-1, coordInf[index][1], :] = [255, 0, 0]

        cv2.imshow('full image', I)
        cv2.waitKey(0) & 0xFF
        cv2.destroyAllWindows()
    
    return manu_pts_dict
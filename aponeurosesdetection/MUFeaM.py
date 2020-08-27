def muscleThickness(I, points1, points2, start, end, calibV, calibH):
    """ Function that calculates muscle thickness  in mm on the 
    interval of columns [start, end] of a US muscle image I.
    
    Args:
        I (array): three-canal image
        points1 (array): dimensions = (n,2)
                        contains the interpolated points of one 
                        aponeurosis in I.
        points2 (array): dimensions = (n,2)
                        contains the interpolated points of the 
                        other aponeurosis in I.
        start (int) : starting column to calculate muscle thickness. Its
                        value should be >= 0 and <I.shape[1]
        end (int): ending column to calculate muscle thickness. Its value
                    should be > start, > 0 and <I.shape[1]
        calibV (float) : vertical calibration factor
        calibH (float) : horizontal calibration factor
    
    Outputs:
            absc (list) : list of abscissa representing distances the beginning
                        of the image I
            mt (list) : muscle thickness in mm at the respecting abscissa
            spl (tuple): spline that interpolates the previous points (x,y)
    """
    if start < 0 or end < 0  or end < start\
        or start >= I.shape[1] or end >= I.shape[1]:
        raise ValueError("start or end values go out of bounds of I's dimensions")

    mt = []
    absc = []
    for col in range(start, end + 1):

        #search if there exist points in each aponeurosis with abscissa col
        search1 = [pt for pt in points1 if pt[1]==col]
        search2 = [pt for pt in points2 if pt[1]==col]

        if search1 and search2: #if there exist a point whith abscissa col 
                                #in each list points1, points2, then calculates MT
            mt.append(abs(search1[0][0]-search2[0][0])*calibV)
            absc.append(col*calibH)
    
    import scipy.interpolate as interpolate
    spl = interpolate.UnivariateSpline(absc, mt, k=5, ext = 0) 
    return absc, mt, spl

def diffSpline(x, spl1, spl2):
    return spl1(x) - spl2(x)

def findIntersections(spl1, listSpl, I, typeI):
    """Function that finds the intersection point between
    two curves respectively defined by the spline spl1
    and the spline spl2

    Args:
        spl1: spline, as output by the scipy function 'univariateSpline'.
                    spl1 accounts for the aponeurosis approximation
        ListSpl (list of splines): list of splines (output by univariateSpline)
                that account for the muscle fascicles approximation

    Outputs:
        listIntersections (list of tuples): intersection points between spl1
        and the different splines of listSpl
    """
    import scipy.optimize as scio
    listIntersections = []
    
    t = 0.01
    if typeI == 'simple':
        x_ini = int(0)
    if typeI == 'panoramic':
        x_ini = int(I.shape[1]/4)
    
    for ind in range(len(listSpl)):
        spl2 = listSpl[ind]
        y0 = int(scio.fsolve(diffSpline, x0 = x_ini, args = (spl1, spl2), xtol = t))
        x0 = int(spl2(y0))
        listIntersections.append((x0,y0))

    return listIntersections

def pennationAngles(spl_a, listS_f, listI, xcalib, ycalib, I = None):
    """
    Function that returns the list of angles between spl_a and
    all curves caracterized by the splines from listS_f

    Args:
        spl_a : spline caracterizing an aponeurosis
        listS_f : list of splines, each one caracterizing
                a muscle fascicle
        listI: list of intersection points between spl_a and 
        splines from listS_f

    Outputs:
        listPA: list of angles in degrees
    """
    import math
    import numpy as np
    list_pa = []

    for index in range(len(listI)):
        I0 = listI[index]
        spl_f = listS_f[index]
        vect = np.arange(I0[1], I0[1] + 150)

        #tangent to fascicle at point I0
        fprime_I0 = (spl_f(I0[1] + 1) - spl_f(I0[1] - 1)) / 2.
        vect_f = fprime_I0 * vect + spl_f(I0[1]) - I0[1] * fprime_I0
        a1 = math.tan( (vect_f[-1] - vect_f[0]) * xcalib / ((vect[-1] - vect[0]) * ycalib) ) * 180 / math.pi

        #tangent to aponeurosis at point I0
        aprime_I0 = (spl_a(I0[1] + 1) - spl_a(I0[1] - 1)) / 2.
        vect_a = aprime_I0 * vect + spl_a(I0[1]) - I0[1] * aprime_I0
        a0 = math.tan( (vect_a[-1] - vect_a[0]) * xcalib / ((vect[-1] - vect[0]) * ycalib) ) * 180 / math.pi

        #Angle computation
        A = abs(a1 - a0)
        list_pa.append(A)

        """
        #Visualize tangents
        for ind in range(vect.shape[0]):
            if vect[ind]>= 0 and vect[ind]<I.shape[1]\
                and vect_f[ind]>=0 and vect_f[ind] < I.shape[0]:
                    I[int(vect_f[ind]), int(vect[ind]), :] = [255,0,255]

            if vect[ind]>= 0 and vect[ind]<I.shape[1]\
                and vect_a[ind]>=0 and vect_a[ind] < I.shape[0]:
                    I[int(vect_a[ind]), int(vect[ind]), :] = [255,0,255]
        """
        
    return list_pa


def fasciclesLength(listS_f, listIu, listId, xcalib, ycalib):
    """
    """
    import numpy as np
    import math
    
    listFL = []
    
    for f in range(len(listS_f)):
        spl = listS_f[f]
        length = 0
        i_u = listIu[f] #intersection point with upper apo
        i_d = listId[f] #intersection point with deep apo
        rangeC = np.arange(min(int(i_d[1]), int(i_u[1])), max(int(i_d[1]), int(i_u[1])) + 1) #columns between these 2 intersection points
        
        for p in range(1, rangeC.shape[0]):
            #compute distance between pixel p and pixel p-1
            deltaX = spl(rangeC[p]) - spl(rangeC[p-1])
            deltaY = rangeC[p]-rangeC[p-1]
            length = length + math.sqrt((deltaX * xcalib) **2 + (deltaY * ycalib) **2)

        listFL.append(length)
        
    return listFL

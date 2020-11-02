"""MUFeaM: MUscle Features Measurements"""
import numpy as np
import math

def pointsCoordinates(typeA, param, interval):
    """
    Creates points coordinates in a given interval,
    according to the chosen type of approximation typeA and its parameters param.

    Inputs:
        typeA (string) : either 'spline' or 'linear'
        param: its type depends on typeA. If typeA = 'spline', then
        param is a spline. If typeA = 'linear', then param is a list of
        2 parameters: [a,b], that caracterises the line equation: row = a*column+b. In
        this case, the line is then converted to a spline and returned by the function
    Outputs:
        an array of the coordinates of the points created from the parameters 
        in the interval (dimensions Nx2)
        param (modified if the type was 'linear', same as input if the type was 'spline')
    """
    if typeA == 'spline':
        y = np.arange(int(interval[0]),int(interval[1]) + 1, 1)
        x = np.int32(param(y))
        coord = np.vstack((x, y)).T  
        
    if typeA == 'linear':
        
        if len(param)<2:
            raise ValueError('Missing parameter for line equation. param should have 2 coefficients.')
        
        y = np.arange(int(interval[0]),int(interval[1]) + 1, 1)
        x = np.int32(y*param[0]+param[1])
        coord = np.vstack((x, y)).T
        
        #transform the equation parameters into a spline
        import scipy.interpolate as interpolate
        param = interpolate.UnivariateSpline(y, x, k=1, ext = 0)
    
    
    return coord, param

def muscleThickness(start, end, calibV, calibH, spl1 = None, spl2 = None, points1 = None, points2 = None):
    """ Function that calculates muscle thickness  in mm on the 
    interval of columns [start, end] of a US muscle image I.
    
    Args:
        spl1 (scipy spline): spline modeling one aponeurosis in I.
        spl2 (scipy spline): spline modeling the other aponeurosis in I.
        points1, points2 (optional): are array of dimension (n,2). 
                they are aponeuroses' points. Default value is None for
                each array (they are calculated from spl1 and spl2 in this case)
        start (int) : starting column to calculate muscle thickness. Its
                        value should be >= 0 and <I.shape[1]
        end (int): ending column to calculate muscle thickness. Its value
                    should be > start, > 0 and <I.shape[1]
        calibV (float) : vertical calibration factor
        calibH (float) : horizontal calibration factor
    
    Outputs:
            absc (list) : list of abscissa representing distances from the beginning
                        of the image I
            mt (list) : muscle thickness in mm at the respecting abscissa
            spl (tuple): spline that interpolates the previous points (x,y)
    """
    start = int(start)
    end = int(end)
    
    mt = []
    absc = []
    
    #generate aponeuroses' points coordinates if necessary
    if (spl1 is None and points1 is None) or (spl2 is None and points2 is None):
        raise ValueError('Missing value. Spline or array of points should be input for each aponeurosis.')
    if points1 is None:
        points1, spl1 = pointsCoordinates('spline', spl1, [start, end])
    if points2 is None:
        points2, spl2 = pointsCoordinates('spline', spl2, [start, end])
    

    for col in range(start, end + 1):
        #search if there exist points in each aponeurosis with abscissa col
        search1 = [pt for pt in points1 if pt[1]==col]
        search2 = [pt for pt in points2 if pt[1]==col]

        if search1 and search2: #if there exist a point whith abscissa col 
                                #in each list points1, points2, then calculates MT
            mt.append(abs(search1[0][0]-search2[0][0])*calibV)
            absc.append(col*calibH)
    
    #interpolation of MT curve
    import scipy.interpolate as interpolate
    spl = interpolate.UnivariateSpline(absc, mt, k=5, ext = 0) 

    return absc, mt, spl



def _diffSpline(x, spl1, spl2):
    """ hidden function that computes the difference of two splines spl1 and spl2 at point x"""
    return spl1(x) - spl2(x)

def findIntersections(spl_inf, spl_sup, listSpl, search_interval, signOfSlope, start = 0):
    """Function that finds the intersection point between
        - spl1 and each spline in the list listSpl
        - spl2 and each spline in the list listSpl
    Conditions to keep considering a fascicle:
        - one intersection with spl_inf and one with spl_sup
        - the orientation of the fascicle must follow the dominant orientation of all fascicles:
            _ if signOfSlope >0:
                - column of intersection with spl_inf > column of intersection with spl_sup
                - line of intersection with spl_inf > line of intersection with spl_sup
            _ if signOfSlope <0:
                - column of intersection with spl_inf < column of intersection with spl_sup
                - line of intersection with spl_inf > line of intersection with spl_sup
        - distance between the two intersection points > 100 pixels (the max calibration factor of all our images is such that 100 pixels = 34mm)
        - intersection points are within the search_interval.
    Otherwise, the fascicle and its intersections with aponeuroses are removed
    from the rest of the analysis
    
    Args:
        spl_inf, spl_sup: splines, as output by the scipy function 'univariateSpline'.
                    they account for the aponeuroses approximations
        ListSpl (list of splines): list of splines (output by univariateSpline)
                that account for the muscle fascicles approximation
        search_interval = list [a,b] = interval in which intersections
                are looked for
        signOfSlope: positive or negative number (integer or float), that
            caracterizes the dominant orientation of fascicles within the muscle
        start (float): parameter for the fsolve function that corresponds to the abscissa
            at which the search for intersection begins.
        
    Outputs:
        listIntersections_i, listIntersections_s (list of tuples): intersection
        points between aponeuroses splines and the different splines of listSpl
        spl_output: new list of splines, where fascicles that do not respect the
        points highlighted above are removed
    """
    import scipy.optimize as scio
    listIntersections_i = []
    listIntersections_s = []
    spl_output = []
    
    t = 0.01
    tol = 0.5
    for ind in range(len(listSpl)):
        spl3 = listSpl[ind]
        a = search_interval[0]
        b = search_interval[1]

        #find columns where interception:
        res1 = scio.fsolve(_diffSpline, x0 = start, args = (spl_inf, spl3), xtol = t, full_output = True)
        res2 = scio.fsolve(_diffSpline, x0 = (b-a)/2, args = (spl_sup, spl3), xtol = t, full_output = True)
        
        #check that solutions have been found and that they are indeed solutions
        if res1[2] == 1 and res2[2] == 1 and abs(res2[1]['fvec'][0])<tol and abs(res1[1]['fvec'][0])<tol:   
            #filter results so that they correspond to the properties of fibres in our images
            if signOfSlope <0:
                if res1[0]<res2[0] and spl3(res1[0])>spl3(res2[0]):
                    #minimal length of fascicle between the 2 intersection points:
                    if ((res1[0]-res2[0])**2 + (spl3(res1[0])-spl3(res2[0]))**2)>100**2:
                        #finally, check if intersections are in the range of the muscle:
                        if res1[0]>a and res1[0]<b and res2[0]>a and res2[0]<b:
                            listIntersections_i.append([int(spl3(res1[0])),int(res1[0])])
                            listIntersections_s.append([int(spl3(res2[0])),int(res2[0])])
                            spl_output.append(spl3)

            elif signOfSlope >0:
                if res1[0]>res2[0] and spl3(res1[0])>spl3(res2[0]):
                    #minimal length of fascicle between the 2 intersection points:
                    if ((res1[0]-res2[0])**2 + (spl3(res1[0])-spl3(res2[0]))**2)>100**2:
                        #finally, check if intersections are in the range of the muscle:
                        if res1[0]>a and res1[0]<b and res2[0]>a and res2[0]<b:
                            listIntersections_i.append([int(spl3(res1[0])),int(res1[0])])
                            listIntersections_s.append([int(spl3(res2[0])),int(res2[0])])
                            spl_output.append(spl3)
         
    return listIntersections_i, listIntersections_s, spl_output


def pennationAngles(spl_a, listS_f, listI, xcalib, ycalib, I = None):
    """
    Function that returns the list of angles between spl_a and
    all curves caracterized by the splines from listS_f

    Args:
        spl_a : spline caracterizing an aponeurosis
        listS_f : list of splines, each one caracterizing a muscle fascicle
        listI: list of intersection points between spl_a and splines from listS_f
        xcalib (float): vertical calibration factor
        ycalib (float): horizontal calibration factor
        I (optional) : three-canal image, to visualize tangents

    Outputs:
        list_pa: list of angles in degrees
    """
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

        if I is not None :
            #Visualize tangents
            for ind in range(vect.shape[0]):
                if vect[ind]>= 0 and vect[ind]<I.shape[1]\
                    and vect_f[ind]>=0 and vect_f[ind] < I.shape[0]:
                    I[int(vect_f[ind]), int(vect[ind]), :] = [255,  0,255]

                if vect[ind]>= 0 and vect[ind]<I.shape[1]\
                    and vect_a[ind]>=0 and vect_a[ind] < I.shape[0]:
                    I[int(vect_a[ind]), int(vect[ind]), :] = [255, 0, 255]
        
    return list_pa


def fasciclesLength(listS_f, listIu, listId, xcalib, ycalib):
    """
    Computes the length of the fascicle between the intersection points with aponeuroses.
    To take into account the potential curvature of our fascicles, the length is the sum
    of the distance between each neighbor pixels within the intersection points.
    
    Inputs:
        listS_f: list of splines that model fascicles
        listIu: list of intersection points of the above spline with superficial aponeurosis
        listId: list of intersection points of the above spline with deep aponeurosis
        xcalib (float): vertical calibration factor
        ycalib (float): horizontal calibration factor

    Outputs:
        list of fascicles length in mm (list has same length as listS_f)
    """
    
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

def locateFasc(intersections, refPoint, ycalib):
    """
    Function that caracterizes each fascicle location by a float value
    corresponding to the horizontal distance between
    refPoint and the intersection of the fascicles with deep aponeurosis
    Ideally, refPoint is the intersection point between the two aponeuroses

    Inputs:
        intersections: list of intersection points of fascicles with deep aponeurosis
        refPoint = (row = on axis 0, column = on axis 1)
        ycalib: horizontal calibration factor
    """
    listLoc = []

    for index in range(len(intersections)):
        pt = intersections[index]
        loc = abs(pt[1] - refPoint[1]) * ycalib
        listLoc.append(loc)

    return listLoc

def curvature(c_points, spline = None):
    """
    # this function has not been tested since most of our fascicles are lines...
    # it should compute the value of the curvature at each point of the input curve
    # I do not know if it works

    inputs:
        spline: function that represents the curve. Default is None
        c_points ( 2D array): 
            if spline is None: array with dimensions Nx2,
            containing N points of the curve. First column is the x
            coordinates (corresponds to the line in the image),
            second column is the y-coordinates ( = columns in image)
            if spline is not None: array of dimensions Nx1,
            which corresponds to the interval of pixels on which 
            curvature must be evaluated

        Outputs
            curva should be an array
    """

    if spline is None:
        import numpy as np
        y_prime = np.gradient(c_points[:,1], edge_order=1,axis=0)
        x_prime = np.gradient(c_points[:,0], edge_order=1,axis=0)
        y_prime2 = np.gradient(y_prime, edge_order=1,axis=0)
        x_prime2 = np.gradient(x_prime, edge_order=1,axis=0)
        curva = (x_prime*y_prime2 - y_prime*x_prime2)/((x_prime**2+y_prime**2)**(3/2))

    elif spline is not None:
        spl_prime = spline.derivative(n=1)
        spl_prime2 = spline.derivative(n=2)
        curva = spl_prime2(c_points[:,1]) / ((1 + spl_prime(c_points[:,1])*spl_prime(c_points[:,1]))**(3/2))

    return curva